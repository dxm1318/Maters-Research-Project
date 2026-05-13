%% Optimization Loop for Sweep Rate (r)

load('PreGenerated_OAE_R.mat'); 
%also loads fs, f1_min, and f1_max

%% Start of Simulation
alpha = 1.2; 
L1_dB = 60; 
L2_dB = 50; 
%fs = 30000;
%f1_min = 700; 
%f1_max = 9000;

sweep_rates = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0];

num_sims = 10;

%vector of distortion products
f_dp_min = [2*f1_min - alpha*f1_min,2*alpha*f1_min-f1_min, 3*f1_min - 2*alpha*f1_min];

% Define the "template" structure with correct field names
template = struct('r', [], 'snr', zeros(1, num_sims), 'latency_ms', zeros(1, num_sims), 'overlap_percent', zeros(1, num_sims));

% Pre-allocate the array using the template
results_r = repmat(template, length(sweep_rates), 3);

for i = 1:length(sweep_rates)
    r = sweep_rates(i);
    L = 1 / (r * log(2)); 
    T = L * log(f1_max / f1_min); 
     
    %predicted latencies of DPs
    dt_dp = -L * log(2-alpha);
    dt_2f2_f1 = -L *log(2*alpha-1);
    dt_3f1_2f2 = -L * log(3-2*alpha); 
    
    dt_vals = [dt_dp,dt_2f2_f1+T,dt_3f1_2f2]*fs; %convert time(s) to samples
    
    for sub = 1:num_sims
        
        %retrieve OAE from dataset
        y = oae_data_r{i,sub};
       
        %identify noise level
        noise_level = 0.5*rms(y);
        %background noise (AWGN)   
        noise = noise_level*randn(size(y));
        %add noise to verhulst output
        y = y + noise;        
        %derive raw impulse response through deconvolution
        [h_raw, t] = VIR_deconv(y,f1_min,L,fs);             
      
        %separate individual ImIRs
        len_IR = 1024;
        pre_IR = 256; 
        hm = synchronized_swept_sine_IR_separation(h_raw, dt_vals, len_IR, pre_IR); 
   
        %separate SL/LL components for each ImIR
        [SL,LL] = separate_SL_LL(hm,fs,r,f_dp_min,T,L,pre_IR);

        for j = 1:length(dt_vals)

            results_r(i,j).r = r;

            metrics_dp = compute_metrics(hm(:,j), SL(:,j), LL(:,j),fs,pre_IR);

            results_r(i,j).snr(sub) = metrics_dp.SNR_SL;
            results_r(i,j).overlap_percent(sub) = metrics_dp.Overlap_percent;
            results_r(i,j).latency_ms(sub) = metrics_dp.latency_ms;          
        end

        fprintf('Simulation %d completed for sweep rate %g\n', sub, r);       
    end
end
fprintf('Simulations are Complete!');
%% Organize Results
organize_results_r(results_r);
        
%% Visualization of Metrics Across Sweep Rates

% Extract sweep rates
r_vals = arrayfun(@(x) x.r, results_r(:,1));
num_rates = length(r_vals);

% Preallocate mean and error arrays
mean_SNR = zeros(num_rates, 3);
se_SNR   = zeros(num_rates, 3);
mean_lat = zeros(num_rates, 3);
se_lat   = zeros(num_rates, 3);
mean_ovlp = zeros(num_rates, 3);
se_ovlp  = zeros(num_rates, 3);

% Compute means and standard errors across simulations
for i = 1:num_rates
    for dp = 1:3
        % Extract data vectors
        snr_data  = results_r(i,dp).snr;
        lat_data  = results_r(i,dp).latency_ms;
        ovlp_data = results_r(i,dp).overlap_percent;
        
        % Number of simulations (n)
        n = length(snr_data);
        
        % SNR Mean and SE
        mean_SNR(i,dp) = mean(snr_data);
        se_SNR(i,dp)   = std(snr_data) / sqrt(n);
        
        % Latency Mean and SE
        mean_lat(i,dp) = mean(lat_data);
        se_lat(i,dp)   = std(lat_data) / sqrt(n);
        
        % Overlap Mean and SE
        mean_ovlp(i,dp) = mean(ovlp_data);
        se_ovlp(i,dp)   = std(ovlp_data) / sqrt(n);
    end
end

%   Plot 1: Sweep Rate vs SNR
% ============================
f1 = figure;
errorbar(r_vals, mean_SNR(:,1), se_SNR(:,1), '-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
errorbar(r_vals, mean_SNR(:,2), se_SNR(:,2), '-s', 'LineWidth', 2, 'MarkerSize', 8);
errorbar(r_vals, mean_SNR(:,3), se_SNR(:,3), '-x', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Sweep Rate r (oct/s)');
ylabel('SNR_{SL} (dB)');
title('SNR of SL Component vs Sweep Rate (Mean \pm SE)');
legend('2f_1 - f_2','2f_2 - f_1','3f_1 - 2f_2', 'Location', 'best');
grid on;
saveas(f1,'SNR_r_noise','png');


%   Plot 2: Sweep Rate vs Latency
% ============================
f2 = figure;
errorbar(r_vals, mean_lat(:,1), se_lat(:,1), '-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
errorbar(r_vals, mean_lat(:,2), se_lat(:,2), '-s', 'LineWidth', 2, 'MarkerSize', 8);
errorbar(r_vals, mean_lat(:,3), se_lat(:,3), '-x', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Sweep Rate r (oct/s)');
ylabel('Separation (ms)');
title('Separation vs Sweep Rate (Mean \pm SE)');
legend('2f_1 - f_2','2f_2 - f_1','3f_1 - 2f_2', 'Location', 'best');
grid on;
saveas(f2,'latency_r_noise','png');

%   Plot 3: Sweep Rate vs Overlap
% ============================
f3 = figure;
errorbar(r_vals, mean_ovlp(:,1), se_ovlp(:,1), '-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
errorbar(r_vals, mean_ovlp(:,2), se_ovlp(:,2), '-s', 'LineWidth', 2, 'MarkerSize', 8);
errorbar(r_vals, mean_ovlp(:,3), se_ovlp(:,3), '-x', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Sweep Rate r (oct/s)');
ylabel('Overlap (%)');
title('SL/LL Overlap vs Sweep Rate (Mean \pm SE)');
legend('2f_1 - f_2','2f_2 - f_1','3f_1 - 2f_2', 'Location', 'best');
grid on;
saveas(f3,'overlap_r_noise','png');

%% Multi-Metric Statistical Analysis (With Interaction)
clf;
close all;

% 1. Load the Data
filename = 'r_sim_Results_noise.csv';
data = readtable(filename);

% Define the output folder
output_folder = 'Analysis_Results';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% 2. Define Metrics and Factors
metrics = {'SNR_SL', 'Latency_ms', 'Overlap__'};
metric_titles = {'SNR SL (dB)', 'Latency (ms)', 'Overlap (%)'};

factors = {data.SweepRate, data.DP_Type};
varnames = {'SweepRate', 'DP_Type'};

fprintf('Starting batch analysis (With Interaction) for %d metrics...\n', length(metrics));

% 3. Loop through each metric
for i = 1:length(metrics)
    current_metric = metrics{i};
    current_title = metric_titles{i};
    response = data.(current_metric);

    fprintf('Processing: %s\n', current_metric);

    % --- Step A: Run ANOVA WITH interaction ---
    [p, tbl, stats] = anovan(response, factors, ...
        'model', 'interaction', ...
        'varnames', varnames, ...
        'display', 'on');

    % --- Step B: Save ANOVA Table to a clean text file ---
    table_filename = fullfile(output_folder, [current_metric, '_ANOVA_Table.txt']);
    fid = fopen(table_filename, 'w');

    fprintf(fid, 'ANOVA Results for %s\n', current_title);
    fprintf(fid, '%s\n\n', repmat('=', 1, 100));

    % Convert ANOVA cell table into strings
    nRows = size(tbl, 1);
    nCols = size(tbl, 2);
    str_tbl = strings(nRows, nCols);

    for r = 1:nRows
        for c = 1:nCols
            val = tbl{r, c};

            if isempty(val)
                str_tbl(r, c) = "";
            elseif ischar(val) || isstring(val)
                str_tbl(r, c) = string(val);
            elseif isnumeric(val)
                if any(isnan(val))
                    str_tbl(r, c) = "";
                else
                    str_tbl(r, c) = sprintf('%.4g', val);
                end
            else
                str_tbl(r, c) = string(val);
            end
        end
    end

    % Calculate automatic column widths
    colWidths = zeros(1, nCols);
    for c = 1:nCols
        colWidths(c) = max(strlength(str_tbl(:, c))) + 4;
    end

    % Print formatted table
    for r = 1:nRows
        for c = 1:nCols
            fprintf(fid, '%-*s', colWidths(c), str_tbl(r, c));
        end
        fprintf(fid, '\n');

        % Separator line after header row
        if r == 1
            for c = 1:nCols
                fprintf(fid, '%s', repmat('-', 1, colWidths(c)));
            end
            fprintf(fid, '\n');
        end
    end

    fclose(fid);

    % --- Step C: Multcompare - Sweep Rate ---
    fig1 = figure;
    multcompare(stats, 'Dimension', 1);
    title(['Main Effect: ', varnames{1}, ' on ', current_title]);
    saveas(fig1, fullfile(output_folder, [current_metric, '_MC_SweepRate.png']));

    % --- Step D: Multcompare - DP Type ---
    fig2 = figure;
    multcompare(stats, 'Dimension', 2);
    title(['Main Effect: ', varnames{2}, ' on ', current_title]);
    saveas(fig2, fullfile(output_folder, [current_metric, '_MC_DPType.png']));

    % --- Step E: Multcompare - Interaction / parameter combinations ---
    fig3 = figure;
    multcompare(stats, 'Dimension', [1 2]);
    title(['Interaction: ', varnames{1}, ' × ', varnames{2}, ' for ', current_title]);
    saveas(fig3, fullfile(output_folder, [current_metric, '_MC_Interaction.png']));
end

fprintf('Analysis complete. ANOVA results saved in: %s\n', output_folder);

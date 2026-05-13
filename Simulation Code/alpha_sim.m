%% Optimization Loop for frequency ratio alpha (f2/f1)

load('PreGenerated_OAE_Alpha.mat');

%% Start of Simulation

%initial parameters
%r = 2.5; %sweep rate oct/s
L1_dB = 60;
L2_dB = 50; 
%fs = 30000; %sample rate
L = 1 / (r*log(2)); %sweep rate coefficient
%f1_min = 700;
%f1_max = 9000;
T = L*log(f1_max/f1_min);
num_sims = 10; %number of simulations

alphas = 1.18:0.01:1.28; %vector of alpha values

%Define template structure and Pre-Allocate results array
template = struct('alpha', [], 'snr', zeros(1, num_sims), 'latency_ms', zeros(1, num_sims), 'overlap_percent', zeros(1, num_sims));
results_a = repmat(template,length(alphas),3);

for i = 1:length(alphas)
    alpha = alphas(i);
    f_dp_min = [2*f1_min-alpha*f1_min,2*alpha*f1_min-f1_min,3*f1_min-2*alpha*f1_min];
    
    %calculate delay time
    dt_dp = -L * log(2-alpha);
    dt_2f2_f1 = -L*log(2*alpha-1);
    dt_3f1_2f2 = -L * log(3-2*alpha);

    dt_vals = [dt_dp, dt_2f2_f1+T,dt_3f1_2f2] * fs; %convert delay time to seconds


    for sub = 1:num_sims
        
        %retrieve OAE data
        y = oae_data_alpha{i,sub};

        %identify noise level
        noise_level = 0.5*rms(y);
        %background noise (AWGN)   
        noise = noise_level*randn(size(y));
        %add noise to verhulst output
        %y = y + noise;       
        %derive IR via deconvolution
        [h_raw,t] = VIR_deconv(y,f1_min,L,fs);

        %separate individual ImIRs
        len_IR = 1024;
        pre_IR = 256;
        hm = synchronized_swept_sine_IR_separation(h_raw,dt_vals,len_IR,pre_IR);

        %separate SL/LL components for each ImIR
        [SL,LL] = separate_SL_LL(hm,fs,r,f_dp_min,T,L,pre_IR);

        for j = 1:length(dt_vals)
            results_a(i,j).alpha = alpha;
            metrics_dp = compute_metrics(hm(:,j),SL(:,j),LL(:,j),fs,pre_IR);

            results_a(i,j).snr(sub) = metrics_dp.SNR_SL;
            results_a(i,j).overlap_percent(sub) = metrics_dp.Overlap_percent;
            results_a(i,j).latency_ms(sub) = metrics_dp.latency_ms;
        end
        fprintf('Simulation %d completed for alpha = %g\n',sub,alpha);
    end
    fprintf('Simulation for alpha = %g is complete. \n',alpha);
end
fprintf('Simulations are Complete!');
%% Oraganize Results
organize_results_a(results_a);
%% Visualization of Metrics Across Alpha (f2/f1)
% Extract alpha values
a_vals = arrayfun(@(x) x.alpha, results_a(:,1));
num_alphas = length(a_vals);

% Preallocate mean and error arrays
mean_SNR  = zeros(num_alphas, 3);
se_SNR    = zeros(num_alphas, 3);
mean_lat  = zeros(num_alphas, 3);
se_lat    = zeros(num_alphas, 3);
mean_ovlp = zeros(num_alphas, 3);
se_ovlp   = zeros(num_alphas, 3);

% Compute means and standard errors across simulations
for i = 1:num_alphas  % Loop through all alpha values
    for dp = 1:3      % Loop through each DP type (2f1-f2, 2f2-f1, 3f1-2f2)
        % Extract data vectors from results_a
        snr_data  = results_a(i,dp).snr;
        lat_data  = results_a(i,dp).latency_ms;
        ovlp_data = results_a(i,dp).overlap_percent;
        
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

%   Plot 1: Alpha vs SNR
% ============================
figure;
errorbar(a_vals, mean_SNR(:,1), se_SNR(:,1), '-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
errorbar(a_vals, mean_SNR(:,2), se_SNR(:,2), '-s', 'LineWidth', 2, 'MarkerSize', 8);
errorbar(a_vals, mean_SNR(:,3), se_SNR(:,3), '-x', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Alpha (\alpha = f_2/f_1)');
ylabel('SNR_{SL} (dB)');
title('SNR of SL Component vs Alpha (Mean \pm SE)');
legend('2f_1 - f_2','2f_2 - f_1','3f_1 - 2f_2', 'Location', 'best');
grid on;

%   Plot 2: Alpha vs Latency
% ============================
figure;
errorbar(a_vals, mean_lat(:,1), se_lat(:,1), '-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
errorbar(a_vals, mean_lat(:,2), se_lat(:,2), '-s', 'LineWidth', 2, 'MarkerSize', 8);
errorbar(a_vals, mean_lat(:,3), se_lat(:,3), '-x', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Alpha (\alpha = f_2/f_1)');
ylabel('Separation (ms)');
title('Separation vs Alpha (Mean \pm SE)');
legend('2f_1 - f_2','2f_2 - f_1','3f_1 - 2f_2', 'Location', 'best');
grid on;

%   Plot 3: Alpha vs Overlap
% ============================
figure;
errorbar(a_vals, mean_ovlp(:,1), se_ovlp(:,1), '-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
errorbar(a_vals, mean_ovlp(:,2), se_ovlp(:,2), '-s', 'LineWidth', 2, 'MarkerSize', 8);
errorbar(a_vals, mean_ovlp(:,3), se_ovlp(:,3), '-x', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Alpha (\alpha = f_2/f_1)');
ylabel('Overlap (%)');
title('SL/LL Overlap vs Alpha (Mean \pm SE)');
legend('2f_1 - f_2','2f_2 - f_1','3f_1 - 2f_2', 'Location', 'best');
grid on;

%% Multi-Metric Statistical Analysis (Alpha)
clf;
close all;

% 1. Load the Data
filename = 'alpha_sim_Results_noise.csv'; 
data = readtable(filename);

% Define the output folder
output_folder = 'Analysis_Results';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% 2. Define Metrics and Factors
metrics = {'SNR_SL', 'Latency_ms','Overlap__'}; 
metric_titles = {'SNR SL (dB)', 'Latency (ms)', 'Overlap (%)'};
factors = {data.alpha, data.DP_Type}; 
varnames = {'alpha', 'DP_Type'};

fprintf('Starting batch analysis for %d metrics...\n', length(metrics));

% 3. Loop through each metric
for i = 1:length(metrics)
    current_metric = metrics{i};
    current_title = metric_titles{i};
    response = data.(current_metric);
    
    fprintf('Processing: %s\n', current_metric);
    
    % --- Step A: Run ANOVA ---
    [p, tbl, stats] = anovan(response, factors, ...
        'model', 'interaction', ...   
        'varnames', varnames, ...
        'display', 'on');

    % --- Step B: Save ANOVA Table ---
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

    % --- Step C: Multcompare - alpha ---
    fig1 = figure;
    multcompare(stats, 'Dimension', 1);
    title(['Main Effect: ', varnames{1}, ' on ', current_title]);
    saveas(fig1, fullfile(output_folder, [current_metric, '_MC_alpha.png']));
  
    % --- Step D: Multcompare - DP Type ---
    fig2 = figure;
    multcompare(stats, 'Dimension', 2);
    title(['Main Effect: ', varnames{2}, ' on ', current_title]);
    saveas(fig2, fullfile(output_folder, [current_metric, '_MC_DPType.png']));
    
    % --- Step E: Multcompare - alpha vs DP Type ---
    fig3 = figure;
    multcompare(stats, 'Dimension', [1 2]);
    title(['Post-hoc: ', varnames{1}, ' vs ', varnames{2}, ' for ', current_title]);
    saveas(fig3, fullfile(output_folder, [current_metric, '_PostHoc_ParamVsParam.png']));   
    
end

fprintf('Analysis complete. Results saved in: %s\n', output_folder);

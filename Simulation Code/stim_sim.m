%% Analysis Loop for Stimulus Strength (Levels) with Latency
% This script processes the pre-generated OAE database and includes Latency metrics.

load('PreGenerated_OAE_Stim.mat'); 
% Variables: oae_database, conditions, labels, alpha, fs, r, f1_min, f1_max

%% Start of Analysis
%alpha = 1.23; 
%r = 2.5; 
%fs = 30000; 
%f1_min = 700; 
%f1_max = 9000;

num_conditions = size(conditions, 1);
num_sims = size(oae_data_stim, 2);

L = 1 / (r * log(2)); 
T = L * log(f1_max / f1_min); 

% Predicted time offsets (samples)
dt_dp      = -L * log(2-alpha);
dt_2f2_f1  = -L * log(2*alpha-1);
dt_3f1_2f2 = -L * log(3-2*alpha); 
dt_vals = [dt_dp, dt_2f2_f1 + T, dt_3f1_2f2] * fs; 

f_dp_min = [2*f1_min - alpha*f1_min, 2*alpha*f1_min - f1_min, 3*f1_min - 2*alpha*f1_min];

% Pre-allocate results with latency_ms field
template = struct('condition_label', "", 'L1', [], 'L2', [], ...
                  'snr', zeros(1, num_sims), 'latency_ms', zeros(1, num_sims), 'overlap_percent', zeros(1, num_sims));
results_L = repmat(template, num_conditions, 3);

for c = 1:num_conditions
    for sub = 1:num_sims

        %retrieve OAE
        y = oae_data_stim{c, sub};
        
        %identify noise level
        noise_level = 0.5*rms(y);
        %background noise (AWGN)   
        noise = noise_level*randn(size(y));
        %add noise to verhulst output
        y = y + noise;  

        %deconvolution
        [h_raw, t] = VIR_deconv(y, f1_min, L, fs);             
      
        %DP Separation
        len_IR = 1024;
        pre_IR = 256; 
        hm = synchronized_swept_sine_IR_separation(h_raw, dt_vals, len_IR, pre_IR); 
   
        %component separation
        [SL, LL] = separate_SL_LL(hm, fs, r, f_dp_min, T, L, pre_IR);

        
        for j = 1:length(dt_vals)
            results_L(c, j).condition_label = labels(c);
            metrics_dp = compute_metrics(hm(:, j), SL(:, j), LL(:, j), fs, pre_IR);

            results_L(c, j).snr(sub) = metrics_dp.SNR_SL;
            results_L(c, j).overlap_percent(sub) = metrics_dp.Overlap_percent;
            results_L(c, j).latency_ms(sub) = metrics_dp.latency_ms; % <--- Latency Metric
        end
        fprintf('Simulation %d for condition %g Complete!\n',sub,c);
    end
    fprintf('Simulations complete for condition %d\n',c);
   
end
%% Organize results
organize_results_L(results_L);
%% Visualization of Metrics Across Stimulus Conditions (Mean ± SE)

% Preallocate mean and error arrays
mean_SNR  = zeros(num_conditions, 3);
se_SNR    = zeros(num_conditions, 3);
mean_lat  = zeros(num_conditions, 3);
se_lat    = zeros(num_conditions, 3);
mean_ovlp = zeros(num_conditions, 3);
se_ovlp   = zeros(num_conditions, 3);

% Compute means and standard errors across simulations
for c = 1:num_conditions
    for dp = 1:3
        % Extract data vectors
        snr_data  = results_L(c,dp).snr;
        lat_data  = results_L(c,dp).latency_ms;
        ovlp_data = results_L(c,dp).overlap_percent;
        
        n = length(snr_data);
        
        % SNR Mean and SE
        mean_SNR(c,dp) = mean(snr_data);
        se_SNR(c,dp)   = std(snr_data) / sqrt(n);
        
        % Latency Mean and SE
        mean_lat(c,dp) = mean(lat_data);
        se_lat(c,dp)   = std(lat_data) / sqrt(n);
        
        % Overlap Mean and SE
        mean_ovlp(c,dp) = mean(ovlp_data);
        se_ovlp(c,dp)   = std(ovlp_data) / sqrt(n);
    end
end

% 1. Use sortrows to sort by L2 (col 2) then L1 (col 1)
% Negative sign '-' denotes descending order
[sorted_conditions, sort_idx] = sortrows(conditions, [-2, -1]);

% 2. Reorder all statistical arrays using the new index
mean_SNR_sorted  = mean_SNR(sort_idx, :);
se_SNR_sorted    = se_SNR(sort_idx, :);
mean_lat_sorted  = mean_lat(sort_idx, :);
se_lat_sorted    = se_lat(sort_idx, :);
mean_ovlp_sorted = mean_ovlp(sort_idx, :);
se_ovlp_sorted   = se_ovlp(sort_idx, :);

% 3. Reorder the labels for the X-axis
sorted_labels = labels(sort_idx);

x_axis = 1:num_conditions;

% --- Now use these 'sorted' variables in your plotting code ---
% (Plotting code remains the same, just ensure you use 'mean_SNR_sorted', etc.)
% Plot 1: Stimulus Condition vs SNR
figure;
errorbar(x_axis, mean_SNR_sorted(:,1), se_SNR_sorted(:,1), '-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
errorbar(x_axis, mean_SNR_sorted(:,2), se_SNR_sorted(:,2), '-s', 'LineWidth', 2, 'MarkerSize', 8);
errorbar(x_axis, mean_SNR_sorted(:,3), se_SNR_sorted(:,3), '-x', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Stimulus Levels (dB SPL)');
ylabel('SNR_{SL} (dB)');
title('SNR vs Stimulus Level');
xticks(x_axis);
xticklabels(sorted_labels);
xtickangle(45);
legend('2f_1 - f_2','2f_2 - f_1','3f_1 - 2f_2');
grid on;

% Plot 2: Stimulus Condition vs Latency
figure;
errorbar(x_axis, mean_lat_sorted(:,1), se_lat_sorted(:,1), '-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
errorbar(x_axis, mean_lat_sorted(:,2), se_lat_sorted(:,2), '-s', 'LineWidth', 2, 'MarkerSize', 8);
errorbar(x_axis, mean_lat_sorted(:,3), se_lat_sorted(:,3), '-x', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Stimulus Levels (dB SPL)');
ylabel('Separation (ms)');
title('Separation vs Stimulus Level');
xticks(x_axis);
xticklabels(sorted_labels);
xtickangle(45);
legend('2f_1 - f_2','2f_2 - f_1','3f_1 - 2f_2');
grid on;

% Plot 3: Stimulus Condition vs Overlap
figure;
errorbar(x_axis, mean_ovlp_sorted(:,1), se_ovlp_sorted(:,1), '-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
errorbar(x_axis, mean_ovlp_sorted(:,2), se_ovlp_sorted(:,2), '-s', 'LineWidth', 2, 'MarkerSize', 8);
errorbar(x_axis, mean_ovlp_sorted(:,3), se_ovlp_sorted(:,3), '-x', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Stimulus Levels (dB SPL)');
ylabel('Overlap (%)');
title('SL/LL Overlap vs Stimulus Level');
xticks(x_axis);
xticklabels(sorted_labels);
xtickangle(45);
legend('2f_1 - f_2','2f_2 - f_1','3f_1 - 2f_2');
grid on;
%% Multi-Metric Statistical Analysis (With Interaction)

% 1. Load the Data
filename = 'L_sim_Results_noisy.csv';
data = readtable(filename);

% Define the output folder
output_folder = 'Analysis_Results_Levels';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% 2. Define Metrics and Factors
metrics = {'SNR_SL', 'Latency_ms', 'Overlap_percent'};
metric_titles = {'SNR SL (dB)', 'Latency (ms)', 'Overlap (%)'};

factors = {data.Condition, data.DP_Type};
varnames = {'Condition', 'DP_Type'};

fprintf('Starting batch analysis for %d metrics...\n', length(metrics));

% 3. Loop through each metric for ANOVA and Multcompare
for i = 1:length(metrics)
    current_metric = metrics{i};
    current_title = metric_titles{i};
    response = data.(current_metric);

    fprintf('Processing: %s\n', current_metric);

    % Run ANOVA WITH interaction
    [p, tbl, stats] = anovan(response, factors, ...
        'model', 'interaction', ...
        'varnames', varnames, ...
        'display', 'off');

    % --- Save ANOVA Table ---
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

    % Multcompare - Condition
    fig1 = figure('Visible', 'on');
    multcompare(stats, 'Dimension', 1);
    title(['Main Effect: ', varnames{1}, ' on ', current_title]);
    saveas(fig1, fullfile(output_folder, [current_metric, '_MC_Condition.png']));

    % Multcompare - DP Type
    fig2 = figure('Visible', 'on');
    multcompare(stats, 'Dimension', 2);
    title(['Main Effect: ', varnames{2}, ' on ', current_title]);
    saveas(fig2, fullfile(output_folder, [current_metric, '_MC_DPType.png']));

    % Multcompare - Interaction
    fig3 = figure('Visible', 'on');
    multcompare(stats, 'Dimension', [1 2]);
    title(['Interaction: ', varnames{1}, ' × ', varnames{2}, ' for ', current_title]);
    saveas(fig3, fullfile(output_folder, [current_metric, '_MC_Interaction.png']));
    
    close([fig1, fig2, fig3]);
end

fprintf('Analysis complete. Results saved in: %s\n', output_folder);


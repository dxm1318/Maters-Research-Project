;%% Pre-Generation Script: generate_oae_data.m
% Run this ONCE to create your "database"
%% amt start
amt_start; %run this command to access the amtoolbox
amt_mex;
%% Base Parameters/constants
fs = 30000; 
f1_min = 700; 
f1_max = 9000;
L1_dB = 60; 
L2_dB = 50;
num_sims = 10;

% Define the ranges you want to pre-calculate
sweep_rates = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0];

% --- Pre-generate for Sweep Rate Variation (r_sim) ---
oae_data_r = cell(length(sweep_rates), num_sims);
alpha = 1.2;

for i = 1:length(sweep_rates)
    r = sweep_rates(i);
    L = 1 / (r * log(2));
    T = L * log(f1_max / f1_min);
    
    for sub = 1:num_sims
        [x_dual,~,~] = SSS_gen(f1_min, f1_max, alpha, T, L1_dB, L2_dB, fs);
        insig = [x_dual(:), x_dual(:)];
        oae_data_r{i, sub} = get_OAE(insig, L1_dB, fs, sub);
        fprintf('OAE %d generated for r = %g\n',sub,r);
    end
    fprintf('Finished OAEs for r = %g\n', r);
end

fprintf('Data Generation is Complete!');
save('PreGenerated_OAE_R.mat', 'oae_data_r', 'sweep_rates', 'fs', 'f1_min', 'f1_max', '-v7.3');
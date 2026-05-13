%% Pre-Generation Script: generate_oae_data_alpha.m
% Run this ONCE to create the OAE database for alpha variation

%% amt_start
amt_start;
amt_mex;
%% Base Parameters/constants
r = 2.5; 
L1_dB = 60;
L2_dB = 50; 
fs = 30000; 
f1_min = 700;
f1_max = 9000;
num_sims = 10; 
%% Data Generation

%array of alphas
alphas = 1.18:0.01:1.28; 

% Pre-allocate cell array to store OAE signals
% Using a cell array because signal lengths might vary slightly with alpha
oae_data_alpha = cell(length(alphas), num_sims);

L = 1 / (r * log(2)); 
T = L * log(f1_max / f1_min);

for i = 1:length(alphas)
    alpha = alphas(i);
    
    for sub = 1:num_sims
        % Generate dual sweep
        [x_dual, ~, ~] = SSS_gen(f1_min, f1_max, alpha, T, L1_dB, L2_dB, fs);
        
        % Prepare input for Verhulst
        insig = [x_dual(:), x_dual(:)];
        
        % Run Verhulst model (The time-consuming part)
        oae_data_alpha{i, sub} = get_OAE(insig, L1_dB, fs, sub);
        fprintf('OAE %d generated for alpha = %g\n',sub,alpha);
    end
    fprintf('Finished OAEs for alpha = %.2f\n', alpha);
end

fprintf('Data Generation is Complete!');

% Save the data to a .mat file
save('PreGenerated_OAE_Alpha.mat', 'oae_data_alpha', 'alphas', 'fs', 'r', 'f1_min', 'f1_max', '-v7.3');

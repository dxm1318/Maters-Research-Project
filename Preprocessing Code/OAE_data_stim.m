%% Pre-Generation Script: generate_oae_stim_levels.m
% Parameters
alpha = 1.2; 
r = 0.5; 
fs = 30000; 
f1_min = 700; 
f1_max = 9000;
num_sims = 10; 
L = 1 / (r * log(2)); 
T = L * log(f1_max / f1_min); 

% Define Level Pairs [L1, L2]
standard_levels = [70, 60; 65, 55; 60, 50; 50, 45; 45, 40];
scissor_levels  = [63, 60; 61, 55; 59, 50; 57, 45; 55, 40];
conditions = [standard_levels; scissor_levels];
num_conditions = size(conditions, 1);

% Pre-allocate Cell Array: [Condition x Simulation]
oae_data_stim = cell(num_conditions, num_sims);

fprintf('Starting OAE generation for %d level conditions...\n', num_conditions);

for c = 1:num_conditions
    L1_curr = conditions(c, 1);
    L2_curr = conditions(c, 2);
    
    fprintf('Condition %d/10: L1=%d, L2=%d\n', c, L1_curr, L2_curr);
    
    for sub = 1:num_sims
        % Generate synchronized swept sines
        [x_dual, ~, ~] = SSS_gen(f1_min, f1_max, alpha, T, L1_curr, L2_curr, fs);
        
        % Run through Verhulst model via get_OAE
        insig = [x_dual(:), x_dual(:)];
        oae_data_stim{c, sub} = get_OAE(insig, L1_curr, fs, sub);
        fprintf('Simulation %d for condition %d complete!',sub,c);
    end
end

fprintf('Data Generation is Complete!');
% Metadata for labeling
labels = ["Std 70/60", "Std 65/55", "Std 60/50", "Std 50/45", "Std 45/40", ...
          "Sci 63/60", "Sci 61/55", "Sci 59/50", "Sci 57/45", "Sci 55/40"];

save('OAE_stim_Levels.mat', 'oae_data_stim', 'conditions', 'labels', ...
     'alpha', 'fs', 'r', 'f1_min', 'f1_max', '-v7.3');

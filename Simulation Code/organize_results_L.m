function flat_results = organize_results_L(results)
% organize_results_L: Organizes stimulus level simulation results into a table
% and exports to CSV.
%
% Define the names of the DPs in the order they appear in the simulation
dp_names = {"2f1-f2", "2f2-f1", "3f1-2f2"};
flatResults = table();

% Iterate through: Conditions (i) and Distortion Products (j)
% size(results, 1) corresponds to num_conditions
% size(results, 2) corresponds to the 3 DP types
for i = 1:size(results, 1)
    for j = 1:size(results, 2)
        
        % Check if the condition label exists to avoid processing empty entries
        if ~isempty(results(i,j).condition_label)
                   
            % Get the number of simulations for this specific condition
            n = length(results(i,j).snr);
            
            % Create a temporary table for this specific DP at this Stimulus Level
            tempTable = table(...
                repmat(results(i,j).condition_label, n, 1), ... % Level Condition Label
                repmat(dp_names{j}, n, 1), ...                  % DP Identifier
                results(i,j).snr', ...                          % SNR (transposed)
                results(i,j).latency_ms', ...                   % Latency (transposed)
                results(i,j).overlap_percent', ...               % Overlap (transposed)
                'VariableNames', {'Condition', 'DP_Type', 'SNR_SL', 'Latency_ms', 'Overlap_percent'});          
           
            % Append to the main table
            flatResults = [flatResults; tempTable];
        end
    end
end

% Save to CSV
writetable(flatResults, 'L_sim_Results_noisy.csv');
fprintf('Export complete! Level results saved\ncl');

end
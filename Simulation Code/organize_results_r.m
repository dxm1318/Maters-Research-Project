
function flat_results = organize_results_r(results)
% Define the names of the DPs in the order they appear in your loop
dp_names = {"2f1-f2","2f2-f1","3f1-2f2"};
flatResults = table();

% Iterate through both dimensions: Sweep Rates (i) and DPs (j)
for i = 1:size(results, 1)
    for j = 1:size(results, 2)
        
        % Only process if the entry isn't empty (in case search_idx was empty)
        if ~isempty(results(i,j).r)
                   
            n = length(results(i,j).snr);
            
            % Create a temporary table for this specific DP at this Sweep Rate
            tempTable = table(...
                repmat(results(i,j).r, n, 1), ...   % Sweep Rate
                repmat(dp_names{j}, n, 1), ...      % DP Identifier
                results(i,j).snr', ...              % SNR (transposed to column)
                results(i,j).latency_ms', ...       % Latency (transposed to column)
                results(i,j).overlap_percent', ...   % Overlap (transposed to column)
                'VariableNames', {'SweepRate', 'DP_Type', 'SNR_SL', 'Latency_ms', 'Overlap_%'});          
           
            % Append to the main table
            flatResults = [flatResults; tempTable];
        end
    end
end

% Save to CSV
writetable(flatResults, 'r_sim_Results_noise.csv');
fprintf('Export complete! Results saved with DP labels.\n');

end
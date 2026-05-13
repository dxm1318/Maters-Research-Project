%% ------------------------------------------------------------------------
function phi = recursive_exp_half(t, Tc)
n = 10; 
phi = zeros(size(t));

% 1. Pre-calculate lambda_n recursively as per footnote 3
gamma = 1;
for k = 1:(n-1)
    gamma = log(gamma + 1);
end
lambda_n = sqrt(gamma);

% 2. Identify positive indices for calculation
idx = (t >= 0);
if any(idx)
    % 3. Apply scaling factor lambda_n
    z = lambda_n * t(idx) / Tc;
    
    % 4. Recursive Gamma calculation
    I = exp(z.^2); % This is Gamma_1
    for k = 1:(n-1) % Loop 9 times to get Gamma_10
        I = exp(I - 1);
    end
    
    phi(idx) = 1 ./ I;
end
end

%% ------------------------------------------------------------------------
% Recursive exponential window
% ------------------------------------------------------------------------
function w = build_recursive_window(t, Tc)
% Half-window for negative times: Tc = 1 ms
Tc_neg = 0.001;

% Compute recursive exponential window φ(T,Tc)
phi_pos = recursive_exp_half(t, Tc);
phi_neg = recursive_exp_half(-t, Tc_neg);

% Combine halves
w = zeros(size(t));
w(t < 0) = phi_neg(t < 0);
w(t >= 0) = phi_pos(t >= 0);

end



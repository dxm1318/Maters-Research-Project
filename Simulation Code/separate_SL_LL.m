function [SL, LL] = separate_SL_LL(hm, fs, r, f_dp_min,T,L,pre_IR)

[N, numDP] = size(hm);
t = (-pre_IR:(N-pre_IR-1))'/fs;

% Novak parameters
a_short = 0.01;
a_full  = 0.05;
b       = -0.8;

% Output arrays
SL = zeros(N, numDP);
LL = zeros(N, numDP);

% Number of windows (scaled with sweep rate)
num_fc = r * 3000;

for k = 1:numDP

        h = hm(:,k);
        
        f_min = f_dp_min(k);
    
        % Time samples along the sweep
        t_fc = linspace(0, T, num_fc);
    
        % Calculate Instantaneous frequencies (Eq. 27)
        fc_vec = f_min * exp(t_fc / L);
    
        % Matrix of ImIRs (Eq. 28)
        Hmat = repmat(h.', num_fc, 1);
    
        % Build window matrices (Eq. 29)
        W_short = zeros(num_fc, N);
        W_full  = zeros(num_fc, N);
    
        for i = 1:num_fc
            fc = fc_vec(i);
    
            Tc_short = (a_short/2) * (fc/1e3)^(b);
            Tc_full  = (a_full/2)  * (fc/1e3)^(b);
    
            W_short(i,:) = build_recursive_window(t, Tc_short);
            W_full(i,:)  = build_recursive_window(t, Tc_full);
        end
    
        % Time Domain windowing
        H_short_mat = Hmat .* W_short;
        H_full_mat  = Hmat .* W_full;
    
        % FFT each row (Eq. 30)
        H_short_fft = fft(H_short_mat, [], 2);
        H_full_fft  = fft(H_full_mat,  [], 2);
    
    
        % Column-wise averaging (Eq. 31)
        SL_complex   = mean(H_short_fft, 1).';
        FULL_complex = mean(H_full_fft,  1).';
    
        % Long-latency = full – short
        LL_complex = FULL_complex - SL_complex;
    
        % Convert to time domain
        SL(:,k) = real(ifft(SL_complex));
        LL(:,k) = real(ifft(LL_complex));
end
end



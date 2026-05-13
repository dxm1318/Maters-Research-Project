function [h_raw,t] = VIR_deconv(y,f1_min,L,fs)
    
    %calculate threhshold spectrum
    X = synchronized_swept_sine_spectra(f1_min,L,fs,length(y));
    Y = (fft(y)./fs);
    Y = Y(:);
    
    %deconvolution
    H = Y./X; 
    H(1) = 0; %prevent infinite/undefined values
    h_raw = ifft(H,'symmetric');
    N = length(h_raw);
    t = (0:N-1)'/fs; %standard time vector
end
function metrics = compute_metrics(hm, SL, LL, fs,pre_IR)
% Computes:
%   1) SNR of the SL component 
%   2) Temporal separation between SL and LL (peak-to-peak latency)
%   3) % overlap between SL and LL -->  mag_LL / (mag_SL+mag_LL)

[N,numDP] = size(hm);
t = (-pre_IR:N-pre_IR-1)'/fs;

for j = 1:numDP
    SLj = SL(:,j);
    LLj = LL(:,j);
    hmj = hm(:,j);
        
    %% Indices for metric Calculations
    [SL_peak,~] = max(abs(SLj)); 
    [LL_peak,~] = max(abs(LLj));
    

    thresh_SL = 0.1*SL_peak; % SL threshold 
    SL_mag = abs(SLj);
    
    thresh_LL = 0.1*LL_peak; % LL threshold
    LL_mag = abs(LLj);
    
    idx_SL = find(SL_mag>thresh_SL,1,'first'); % approx. SL start
    idx_LL = find(LL_mag>thresh_LL,1,'first'); %approx. LL start
    

    idx1 = idx_SL;
    idx2 = find(t==0.002); %2ms endpoint
    
    %% Approximate Latency between SL and LL components
    latency = t(idx_LL)-t(idx_SL);
    latency_ms = latency * 1000;
    %% SNR SL

    s_rms = rms(SL_mag(idx1:idx2));
    n_rms = rms(hmj(end/2:end)); 
 
    SNR_SL = 20*log10(s_rms/n_rms);
   
    %% Overlap: How much do the the SL and LL signals overlap in time? 

    SL_dist = idx2-idx_SL; %sample distance SL
    LL_dist = idx2-idx_LL; %sample distance LL
    overlap = LL_dist/SL_dist;
    overlap_percent = overlap * 100;
  
    %% Store results
    metrics(j).SNR_SL = SNR_SL;
    metrics(j).Overlap_percent = overlap_percent;
    metrics(j).latency_ms = latency_ms;
    
end
end

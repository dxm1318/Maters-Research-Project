function [x_dual,x1,x2] = SSS_gen(f1_min,f1_max,alpha,T,L1_dB,L2_dB,fs)
    
    %use frequency ratio (alpha) to create f2 frequency range
    f2_min = alpha*f1_min; 
    f2_max = alpha*f1_max;
             
    %individual sweeps
    [x1,~,~] = synchronized_swept_sine(f1_min,f1_max,T,fs); %f1 sweep
    [x2,~,~] = synchronized_swept_sine(f2_min,f2_max,T,fs); %f2 sweep
    
    %generate dual sweep
    rel_f2_dB = L2_dB - L1_dB;
    x_dual = x1 + db2mag(rel_f2_dB) * x2; %x2 is scaled in amplitude

end
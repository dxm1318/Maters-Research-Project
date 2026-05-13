function y = get_OAE(insig,L1_dB,fs, sub)

    %generate virtual OAE using verhulst cochlear model
    [~,~,OAE,~] = verhulst2012(insig,fs,'all', [L1_dB L1_dB],'irr',[1 0], 'subject', sub);
    
    %remove stimulus
    y = OAE(:,1)- OAE(:,2);
    %y = OAE(:,2);
    

end
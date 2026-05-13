function hm = synchronized_swept_sine_IR_separation(h,dt,len_IR,pre_IR)
% SYNCHRONIZED_SWEPT_SINE_IR_SEPARATION - separation of impulse responses
% Separate impulse responses hm, corresponding to higher orders, from the
% impulse response h.
%
%    h = sum_m h_m(t + dt_m)
%
% Syntax:  hm = synchronized_swept_sine_IR_separation(h,dt,len_IR,pre_IR)
%
% Inputs:
%    h      - (complete) impulse response (column vector)
%    dt     - time delays at which the higher harmonics appears
%    len_IR - length of the impulse responses hm [samples]
%    pre_IR - number of samples before the IR
%
% Outputs:
%    hm - separated impulse responses (N-column vector, N=length(dt))
%

% Author: Antonin NOVAK
% Laboratoire d'Acoustique de l'Universite du Mans (LAUM) - UMR CNRS n.6613
% email address: antonin.novak@univ-lemans.fr  
% Website: https://ant-novak.com
% July 2016; Last revision: 29-Oct-2018

%------------- BEGIN CODE --------------

% number of IRs to be separated
N = length(dt);  

% rounded positions [samples]
dt_ = round(dt);

% non-integer difference
dt_rem = dt - dt_;

% circular periodisation of IR
h_pos = [h; h(1:len_IR)];

% frequency axis definition (0 .. 2pi)
axe_w = linspace(0,2*pi,len_IR+1)'; 
axe_w(end) = [];

% memory localisation for IRs
hm = zeros(len_IR,N);

% last sample poition
%st0 = length(h_pos);

for n = 1:N
    % start position of n-th IR
    st = dt_(n) - pre_IR;
    % end position of n-th IR
    ed = st + len_IR;    
    % separation of IRs
    hm(1:ed-st,n) = h_pos(st+1:ed);
    % Non-integer sample delay correction
    Hx = fft(hm(:,n)) .* exp(-1i*dt_rem(n)*axe_w);
    hm(:,n) = ifft(Hx,'symmetric');
    
    % last sample poition for the next IR
    %st0 = st - 1;
end

%------------- END CODE --------------
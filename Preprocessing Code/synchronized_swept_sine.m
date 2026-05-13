function [x,t,L] = synchronized_swept_sine(f1,f2,T,fs)
%SYNCHRONIZED_SWEPT_SINE Synchronized Swept-Sine generator.
% Generates samples of a synchronized swept-sine signal, generated using 
% the following equation [1]
%        x = sin(2*pi*f1*L*(exp(t/L))),
%   where t is the time support and
%        L = T/log(f2/f1)
%   is the rate of frequency increase.
%
% Syntax:  [x,t,L] = SYNCHRONIZED_SWEPT_SINE(f1,f2,T,fs)
%
% Inputs:
%    f1 - instantaneous frequency at time 0 [Hz]
%    f2 - instantaneous frequency at time T [Hz]
%    T  - duration of the swept-sine [seconds]
%    fs - sampling frequency [Hz]
%
% Outputs:
%    s - Synchronized Swept-Sine signal (column vector)
%    t - time support (column vector)
%    L - the rate of frequency increase
%
%   References:
%   [1] Novak, A., et al. "Synchronized Swept-Sine: Theory, Application, and Implementation." Journal of the Audio Engineering Society 63.10 (2015): 786-798.
%   [2] Farina, A. "Simultaneous measurement of impulse response and distortion with a swept-sine technique." Audio Engineering Society Convention 108. Audio Engineering Society, 2000.
%   [3] Novak, A., et al. "Nonparametric identification of nonlinear systems in series." IEEE Transactions on Instrumentation and Measurement 63.8 (2014): 2044-2051.

% Author: Antonin NOVAK
% Laboratoire d'Acoustique de l'Universite du Mans (LAUM) - UMR CNRS n.6613
% email address: antonin.novak@univ-lemans.fr  
% Website: https://ant-novak.com
% July 2016; Last revision: 29-Oct-2018

%------------- BEGIN CODE --------------

L = T/log(f2/f1);
t = (0:1/fs:T-1/fs)';
x = sin(2*pi*f1*L*(exp(t/L)));

%------------- END OF CODE -------------




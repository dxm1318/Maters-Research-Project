function [X,f] = synchronized_swept_sine_spectra(f1,L,fs,Npts)
%SYNCHRONIZED_SWEPT_SINE_SPECTRA Synchronized Swept-Sine spectra.
% Generates the spectra of the synchronized swept-sine signal (for positive 
% frequencies) using the following equation [1]
% X = 1/2*sqrt(L./f_ax).*exp(1i*2*pi*f_ax*L.*(1 - log(f_ax/f1)) - 1i*pi/4),
% where f_ax is the freqeuncy axis and L is the rate of frequency increase.
%
% Syntax:  [X,f_ax] = synchronized_swept_sine_spectra(f1,L,fs,Npts)
%
% Inputs:
%    f1   - instantaneous frequency at time 0 [Hz]
%    f2   - instantaneous frequency at time T [Hz]
%    T    - duration of the swept-sine [seconds]
%    fs   - sampling frequency [Hz]
%    Npts - number of points (length of the spectra)
%
% Outputs:
%    X    - Synchronized Swept-Sine spectra (column vector)
%    f    - frequency axes

%
% References:
%  [1] Novak, A., et al. "Synchronized Swept-Sine: Theory, Application, and Implementation." Journal of the Audio Engineering Society 63.10 (2015): 786-798.

% Author: Antonin NOVAK
%  Laboratoire d'Acoustique de l'Universite du Mans (LAUM) - UMR CNRS n.6613
%  email address: antonin.novak@univ-lemans.fr  
%  Website: https://ant-novak.com
% July 2016; Last revision: 29-Oct-2018

%------------- BEGIN CODE --------------

% frequency axis
f = linspace(0,fs,Npts+1)'; 
f(end) = [];

% spectra of the synchronized swept-sine signal [1, Eq.(42)]
X = 1/2*sqrt(L./f).*exp(1i*2*pi*f*L.*(1 - log(f/f1)) - 1i*pi/4); 
X(1) = Inf; % protect from devision by zero

%------------- END OF CODE -------------




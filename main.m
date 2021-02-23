lambda = 1;  % in meters
L = 0.47*lambda;  % length of the wire
a = 0.005*lambda;  % radius of the wire
nSegments = 11;  % number of segments on the wire
excitedSeg = 6;  % segment where excitation is given
freq = 50;  
k = 2*pi/lambda;  % wave number :: (2*pi)/lambda
V = 1;  % voltage applied on the segment

[coeff, Zin] = wireDipoleAntenna(L, a, nSegments, excitedSeg, freq, k, V);

% plot Xin, Rin vs number of basis functions N
% plot abs(current distribution) vs (z/lambda)
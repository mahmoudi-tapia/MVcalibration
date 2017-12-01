%% Last update: 03/03/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% This function gives the mu** for the MV t-distribution

function [ mu, T ] = muDstar2( x, BetaHAT, D, Achol, H, input, r)

p = size(x,2);
Hx = designMat(x, p+1);

% Get matrix T made of vectors t
Xphi1 = bsxfun(@times,input,sqrt(r'));
Xphi2 = bsxfun(@times,x,sqrt(r'));
Q = bsxfun(@plus,dot(Xphi2,Xphi2,2)',dot(Xphi1,Xphi1,2)) - 2*(Xphi1*Xphi2');
T = exp(-Q);

temp1 = Achol'\(D - H * BetaHAT);
temp2 = Achol'\T;

mu = Hx*BetaHAT + temp2' * temp1;
end
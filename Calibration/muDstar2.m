%% Last update: 03/03/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% This function gives the mu** for the MV t-distribution

function [ mu, T ] = muDstar2( xE, theta, BetaHAT, D, A, H, input, r)

N =size(xE,1);
xE = [xE,repmat(theta',N,1)];
p = size(xE,2);
HxE = designMat(xE, p+1);

% Get matrix T made of vectors t
Xphi1 = bsxfun(@times,input,sqrt(r'));
Xphi2 = bsxfun(@times,xE,sqrt(r'));
Q = bsxfun(@plus,dot(Xphi2,Xphi2,2)',dot(Xphi1,Xphi1,2)) - 2*(Xphi1*Xphi2');
T = exp(-Q);

U = chol(A);
temp1 = U'\(D - H * BetaHAT);
temp2 = U'\T;

mu = HxE*BetaHAT + temp2' * temp1;
end
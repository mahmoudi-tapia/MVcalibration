%% Last update: 04/19/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% This function gives the Sigma_hat_GLS

function [ SigmaHAT_GLS ] = SigmaHAT2( Achol, H, D, B_GLS)

n = size(Achol,1);
m = size(H,2);

temp1 = (D-H*B_GLS);
W1 = Achol'\temp1;
SigmaHAT_GLS = 1/(n-m) * (W1' * W1);

end
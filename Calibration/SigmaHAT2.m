%% Last update: 04/19/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% This function gives the Sigma_hat_GLS

function [ SigmaHAT_GLS ] = SigmaHAT2( A, H, D, B_GLS)

n = size(A,1);
m = size(H,2);
U = chol(A);
W1 = (D-H*B_GLS);
temp1 = U'\W1;
SigmaHAT_GLS = 1/(n-m) * (temp1'*temp1);

end
%% Last update: 03/02/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% This function gives the Beta_hat

function [ BetaHAT_GLS ] = BetaHAT2( A, H, D)

Q1 = chol(A);
W1 = Q1'\H;     % W_1 = Q_1^{-T} * H
Z1 = W1'*W1;   % Z_1 = H^T * A^(-1) * H
Z2 = Z1\H';     % Z_2 = (H^T * A^(-1) * H)^(-1) * H^T

Z3 = Q1'\Z2';
Z4 = Q1'\D;

BetaHAT_GLS = Z3'* Z4;
end
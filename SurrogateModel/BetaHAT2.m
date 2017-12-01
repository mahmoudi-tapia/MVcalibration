%% Last update: 09/01/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% This function gives the Beta_hat

function [ BetaHAT_GLS ] = BetaHAT2( Achol, H, D)

W1 = Achol'\H;     % W_1 = Q_1^{-T} * H
Z1 = W1'*W1;   % Z_1 = H^T * A^(-1) * H
Z2 = Z1\H';     % Z_2 = (H^T * A^(-1) * H)^(-1) * H^T

Z3 = Achol'\Z2';
Z4 = Achol'\D;

BetaHAT_GLS = Z3'* Z4;
end
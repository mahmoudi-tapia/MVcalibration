%% Last update: 09/01/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% This function computes the posterior probability of the vector r 
% which is the vector of positive roughness parameters

function [ postProb ] = log_pi_r_D2( r, D, A, Achol, H, q )

n = size(A,1);
m = size(H,2);

Ainv = inv(A);

W1 = Achol'\H;     % W_1 = Q_1^{-T} H
Z1 = W1'*W1;   % Z_1 H^T * A^(-1) * H


Q2 = chol(Z1);
Z2 = (A\H)';
W2 = Q2' \ Z2;
Z3 = W2'*W2;


% G = temp1 - temp1 * H * inv(Z) * H' * temp1;
G = Ainv - Z3;

% postProb = pi_r(r) * det(A)^(-q/2) * det(temp2)^(-q/2) * det(D'*G*D)^(-(n-m)/2);
%**** Gustavo: changed function pi_r for log_pi_r
postProb = log_pi_r(r) + ...
    (-q/2)*2*sum(log(diag(Achol))) + ...
    (-q/2)*2*sum(log(diag(Q2))) +...
    (-(n-m)/2)*log(det(D'*G*D));

end


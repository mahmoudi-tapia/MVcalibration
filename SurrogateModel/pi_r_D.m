%% Last update: 02/22/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% This function computes the posterior probability of the vector r 
% which is the vector of positive roughness parameters

function [ postProb ] = pi_r_D( r, D, A, H, q )

n = size(A,1);
m = size(H,2);

temp1 = inv(A);
% temp2 = H' * temp1 * H; % without Choleski decomposition
Q = chol(A);
W = Q'\H;
temp2 = W'*W;
G = temp1 - temp1 * H * inv(temp2) * H' * temp1;



postProb = pi_r(r) * det(A)^(-q/2) * det(temp2)^(-q/2) * det(D'*G*D)^(-(n-m)/2);

end


%% Last update: 02/28/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% This function gives the correlation matrix A


function [ A ] = corrMat2( r, input )

% Sanity check
if length(r) ~= size(input, 2)
    disp('BAD INPUT! -- CHECK THE CODE')
    return
end

% Get the size of the input vectors
n = size(input, 1);

% Declare the matrix full of zeros
% A = zeros(n,n);

Xphi = bsxfun(@times,input,sqrt(r'));
Q = bsxfun(@plus,dot(Xphi,Xphi,2)',dot(Xphi,Xphi,2)) - 2*(Xphi*Xphi');
Q(eye(n)==1) = 0;
A = exp(-Q);
end
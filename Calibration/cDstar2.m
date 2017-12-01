%% Last update: 04/19/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% This function gives the c** 

% h(x)
% c(x1, x2)
% t(x)
% Z1

function [ corr ] = cDstar2( xE, theta, T, A, H, r)

N = size(xE,1);
xE_extended = [xE repmat(theta', N, 1)];
p = size(xE_extended,2);
HxE = designMat(xE_extended, p+1);
U = chol(A);

temp1 = corrMat2( r, xE_extended );


W1 = U'\T;
temp2 = W1'*W1;


W2 = U'\H;
W3 = W2'*W1;
temp3 = HxE' - W3;

U2 = chol(W2'*W2);
W4 = U2'\temp3;
temp4 = W4'*W4;

corr = temp1 - temp2 + temp4;

end


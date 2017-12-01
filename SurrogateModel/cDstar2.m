%% Last update: 09/01/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% This function gives the c** 

% h(x)
% c(x1, x2)
% t(x)
% Z1

function [ corr ] = cDstar2( x, T, Achol, H, r)

p = size(x,2);
Hx = designMat(x, p+1);

temp1 = corrMat2( r, x );

W1 = Achol'\T;
temp2 = W1'*W1;

W2 = Achol'\H;
W3 = W2'*W1;
temp3 = Hx' - W3;

U2 = chol(W2'*W2);
W4 = U2'\temp3;
temp4 = W4'*W4;

corr = temp1 - temp2 + temp4;

end


%% Last update: 04/19/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% This function gives the c** 

% h(x)
% c(x1, x2)
% t(x)
% Z1

function [ C_em_PE, C_delta_PE  ] = corrPE( xP, xE, thetaST, T_xP, T_xE, A, H, r_hat, rDeltaST)

N = size(xE,1);
S = size(xP,1);

xE_extended = [xE repmat(thetaST', N, 1)];
xP_extended = [xP repmat(thetaST', S, 1)];

p = size(xE_extended,2);

HxP = designMat(xP_extended, p+1);
HxE = designMat(xE_extended, p+1);


% Get vectors for binary operations for C_delta_PE
XphiP = bsxfun(@times,xP,sqrt(rDeltaST'));
XphiE = bsxfun(@times,xE,sqrt(rDeltaST'));

% Compute C_delta_PE
Q = bsxfun(@plus,dot(XphiE,XphiE,2)',dot(XphiP,XphiP,2)) - 2*(XphiP*XphiE');
C_delta_PE = exp(-Q);


% Get vectors for binary operations for C_em_PE
XphiP = bsxfun(@times,xP_extended,sqrt(r_hat'));
XphiE = bsxfun(@times,xE_extended,sqrt(r_hat'));


% Compute the terms for C_em_PE
Q = bsxfun(@plus,dot(XphiE,XphiE,2)',dot(XphiP,XphiP,2)) - 2*(XphiP*XphiE');
temp1 = exp(-Q);

U = chol(A);

W1E = U'\T_xE;
W1P = U'\T_xP;
temp2 = W1P'*W1E;


W2 = U'\H;

W3E = W2'*W1E;
W3P = W2'*W1P;

temp3E = HxE' - W3E;
temp3P = HxP' - W3P;

U2 = chol(W2'*W2);

W4E = U2'\temp3E;
W4P = U2'\temp3P;

temp4 = W4P'*W4E;

C_em_PE = temp1 - temp2 + temp4;

end


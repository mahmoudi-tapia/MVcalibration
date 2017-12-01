%% Last update: 07/13/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Calibrator

% This function computes the logarithm of the posterior probability for all the parameters
% xE: Input matrix of physical observations (N \times \kappa)
% yE: Output matrix of physical observations (N \times q)

% xS: Input matrix of numerical observations (n \times p)
% D: Output matrix of numerical observations (n \times q)


function [ logP ] = post_PHI2(  theta, r_delta, psi_delta, psi_eps,...
                                xE, yE,...
                                r_hat, xS, D, A, H, BetaHAT)

% Get the model parameters
% Refer to the code implenetation document for calibration
% kappa: number of control variables (laser power and speed)
% N: number of experimental observations
% q: number of outputs
[N, kappa] = size(xE);
[N, q] = size (yE);

% Stack the input and output matrices to make vectors
% xEstacked = reshape(xE', [N*kappa,1]);
yEstacked = reshape(yE', [N*q,1]);

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Compute the mean vector mu^{**} for the MVN which is of size (Nq)x1
% muDS = zeros(N*q,1);
% for i=1:N
%     % Create a single input using calibration data and the parameters
%     x = [xE(i,:), theta']';
%     
%     % Compute assign the right portion of the mean vector
%     muDS(1+(i-1)*q: i*q) = muDstar( x, BetaHAT, D, A, H, xS, r_hat);
% end

[muDS, T] = muDstar2( xE, theta, BetaHAT, D, A, H, xS, r_hat);
muDS = reshape(muDS',N*q,1);
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

%% Compute the covariance matrix for the MVN which is of size (Nq)x(Nq)
% Combine the input matrix of experimental observations with calibration
% parameters
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% xE_extended = [xE repmat(theta', N, 1)];
% 
% % Covariance matrix Part I (emulator)
% C_em=zeros(N,N);
% for i=1:N
%     
%     % Get the first vector
%     xi = xE_extended(i, :)';
%     
%     for j=i:N
%         % Get the second vector
%         xj = xE_extended(j, :)';
%         
%         % Generate (i,j)th element of the matrix
%         temp1 = cDstar(xi, xj, A, H, xS, r_hat);
%         C_em(i,j) = temp1;
%         C_em(j,i) = temp1;
%     end
% end

C_em = cDstar2( xE, theta, T, A, H, r_hat);

SigmaGLS = SigmaHAT2( A, H, D, BetaHAT);
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
covI = kron(C_em, SigmaGLS);

% Covariance matrix Part II (discrepancy)
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C_delta = corrMat2( r_delta, xE); % Note that the discrepancy doesn't depend on calibration parameters hence only xE is input not xE_extended
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
covII = kron(C_delta, diag(psi_delta));

% Covariance matrix Part III (error)
covIII = kron(eye(N), diag(psi_eps));

covMVN = covI + covII + covIII;

% Compute the likelihood
logLike = logmvnpdf(yEstacked', muDS', covMVN);

% Compute the prior probability
logPrior = pi_PHI( theta, r_delta, psi_delta, psi_eps );

% Compute the log posterior prob
logP = logPrior + logLike;
end
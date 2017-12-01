%% Last update: 06/13/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Calibrator

% This function computes the logarithm of the prior probability for all the parameters
% Distribution functions that are used are:
% Uniform for theta
% Log-logistic for r_delta
% Inverse gamma for psi_delta
% Inverse Wishart for psi_eps


function [ logP ] = pi_PHI( theta, r_delta, psi_delta, psi_eps )

% Get the model parameters
% Refer to the code implenetation document for calibration
% kappa: number of control variables (laser power and speed)
% p: number of inputs (control variables + calibration parameters)
% q: number of outputs

q = length(psi_eps);
kappa = length(r_delta);
p = kappa + length(theta);

% In this code we compute four prior probabilities independently

% ----------------------------
% ----------   P1
% ----------------------------
% PARAMETERS for theta
LB=zeros(3,1);    % lower bound for each calib. parameter
UB=ones(3,1);    % upper bound for each calib. parameter

P1=zeros(p-kappa,1);  % define a vector for each probability
for i=1:p-kappa
    P1(i) = unifpdf(theta(i),LB(i), UB(i));
end

P1 = sum(log(P1));  % Joint prob for calib. parameters

% ----------------------------
% ----------   P4
% ----------------------------
% PARAMETERS for psi_eps
IG_A_eps = 2;   % shape parameter
IG_B_eps = 1;   % scale parameter
% Prior probability for psi_eps using built-in function
P4 = sum(log(gampdf(psi_eps,IG_A_eps,IG_B_eps)));

% ----------------------------
% ----------   P2
% ----------------------------
% Use the log-logistic prior for the roughness parameters
P2 = pi_r_delta (r_delta);

% ----------------------------
% ----------   P3
% ----------------------------
% PARAMETERS for psi_delta
IG_A_delta = 2;   % shape parameter
IG_B_delta = 1;   % scale parameter
% Prior probability for psi_eps using built-in function
P3 = sum(log(gampdf(psi_delta,IG_A_delta,IG_B_delta)));






% -------------------- Code 2 (NOT WORKING)
% % PARAMETERS for Sigma_delta
% IW_d = q + 4;  % degree of freedom
% IW_S = [    .65   .5  .5      % positive def. scale matrix
%             .5  .75 .5
%             .5  .5  .75];
% IW_S = cov(yE);
% IW_S = eye(q);
% % Prior probability for Sigma_delta
% [ dummy, P3 ] = invwishpdf(Sigma_delta,IW_S,IW_d);

% -------------------- Code 1 (UNNECESSARY)
% % This code is for a 3-by-1 output, hence Sigma is 3-by-3
% % Sanity check
% if length(psi_delta) ~= (q*(q+1)/2)
%     disp('BAD INPUT! -- CHECK THE CODE')
%     return
% end      
%   
% % Make a symmetric matrix from the input vector psi_delta
% IW = squareform(psi_delta(q+1:end));    % make a square matrix
% IW = IW + diag(psi_delta(1:q));         % set the diagonal





% ----------------------------
% Joint prior probabilities
logP = P1 + P2 + P3 + P4;

end
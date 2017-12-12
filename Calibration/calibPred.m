%% Last update: 11/29/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% This script generates the emulator
    

% SAMPLE FUNCTION RECALL 
% [ mvCalibPred ] = calibPred( r_hat, 'surrogate_model_data.csv', 'calibration_data.csv', 5, 3, 2, mvCalibRes.thetaMode , mvCalibRes.rDeltaHAT, mvCalibRes.psiDeltaHAT, mvCalibRes.psiEpsHAT, 'calib_pred_data.csv');

function [ mvCalibPred ] = calibPred( r_hat, emulTrainData, calibTrainData, P, q, t, thetaMode, rDeltaST, psiDeltaST, psiEpsST, predDataCSV)


% Description of the input arguments

% r_hat             Vector of size (p+kappa)x1 which is the estimate of the surrogate model paramaters
% emulTrainData     numerical experiments (simulation) data filename string - CSV with header
% calibTrainData     physical experiments data filename string - CSV with header
% p                 number of input parameters (control variables + calibration parameters)
% q                 number of output parameters
% t                 number of control inputs
% thetaMode         Mode of the posterior distribution for the calibration parameters
% rDeltaST          Estimate of the roughness parameters of the model outputs
% psiDeltaST        Estimate of the variance of the model outputs
% psiEpsST          Estimate of the noise variance
% predDataCSV       data filename string, Inputs for the purpose of prediction using the emulator after calibrating parameters




%-------------------
% Read the csv files
%-------------------

% Prediction input
predData = csvread(predDataCSV, 1, 0);  % Ignore the file header
S = size(predData,1);

% the calibration data
fullCalibData = csvread(calibTrainData, 1, 0);  % Ignore the file header
exprData = fullCalibData;
N = size(exprData,1);

% the simulation data
simData = csvread(emulTrainData, 1, 0);
xS_RAW = simData(:,1:P);    % Collect the input
% Normalize the input -- between 0 and 1
xS = (xS_RAW-repmat(min(xS_RAW),length(xS_RAW),1))...
            ./repmat([max(xS_RAW)-min(xS_RAW)],length(xS_RAW),1);

D_RAW = simData(:,P+1:P+q);    % Collect the output matrix
% Normalize output -- mean 0
D = (D_RAW-repmat(mean(D_RAW),length(D_RAW),1))...
            ./repmat(std(D_RAW),length(D_RAW),1);   


% Normalize the data  

xP_RAW = predData(:,1:t);    %input
% Normalize the input -- between 0 and 1
temp = xS_RAW(:,1:t);   % Normalize according to the simulation data
xP = (xP_RAW-repmat(min(temp),length(xP_RAW),1))...
            ./repmat([max(temp)-min(temp)],length(xP_RAW),1);

        

        
        
% Collect the input
xE_RAW = exprData(:,1:t);    
% Normalize the input -- between 0 and 1
temp = xS_RAW(:,1:t);   % Normalize according to the simulation data
xE = (xE_RAW-repmat(min(temp),length(xE_RAW),1))...
            ./repmat([max(temp)-min(temp)],length(xE_RAW),1);
        
% Collect the output matrix
yE_RAW = exprData(:,t+1:t+q);    
% Normalize output -- mean 0
yE = (yE_RAW-repmat(mean(D_RAW),length(yE_RAW),1))...
            ./repmat(std(D_RAW),length(yE_RAW),1);    

        
        
%% Prediction
% Generate the design matrix H, the correlation matrix A, and the GLS regression coefficients BetaHAT
% We need the results of the emulation (r_hat) for this step
H = designMat (xS, P+1);
A = corrMat2(r_hat, xS);        
BetaHAT_GLS = BetaHAT2( A, H, D);

thetaST = thetaMode;

% Computing terms neeeded for predicting E[yP|yE]
[muDSxP, T_xP] = muDstar2( xP, thetaST, BetaHAT_GLS, D, A, H, xS, r_hat);
[muDSxE, T_xE] = muDstar2( xE, thetaST, BetaHAT_GLS, D, A, H, xS, r_hat);

% Make the vectors stacked
muDSxP = reshape(muDSxP', [S*q,1]);
muDSxE = reshape(muDSxE', [N*q,1]);
yEstacked = reshape(yE', [N*q,1]);

% Compute SigmaEE
C_em = cDstar2( xE, thetaST, T_xE, A, H, r_hat);
SigmaGLS = SigmaHAT2( A, H, D, BetaHAT_GLS);
covI = kron(C_em, SigmaGLS);    % Term I for SigmaEE

C_delta = corrMat2( rDeltaST, xE); 
covII = kron(C_delta, diag(psiDeltaST));    % Term II for SigmaEE

covIII = kron(eye(N), diag(psiEpsST));  % Term III for SigmaEE

SigmaEE = covI + covII + covIII;

% Compute SigmaPE
[ C_em_PE, C_delta_PE  ] = corrPE( xP, xE, thetaST, T_xP, T_xE, A, H, r_hat, rDeltaST);
SigmaPE = kron(C_em_PE,SigmaGLS) + kron(C_delta_PE,diag(psiDeltaST));

% Compute SigmaPP
[ C_em_PP, C_delta_PP  ] = corrPE( xP, xP, thetaST, T_xP, T_xP, A, H, r_hat, rDeltaST);
SigmaPP = kron(C_em_PP,SigmaGLS) + kron(C_delta_PP,diag(psiDeltaST)) ...
    + kron(eye(S),diag(psiEpsST));



% Kriging for prediction
temp = SigmaEE\(yEstacked - muDSxE);
E_yP = muDSxP + SigmaPE * temp;
E_yP = reshape(E_yP,[q,S])';

% Compute the variance
temp2 = SigmaEE\SigmaPE';
V_yP = SigmaPP - SigmaPE * temp2;

E_yP_actual = E_yP .* repmat(std(D_RAW),S,1) + repmat(mean(D_RAW),S,1);
V_yP_actual = reshape(diag(V_yP),q,S)' .* repmat(std(D_RAW),S,1);

%% Save the resutls in a struct
f1 = 'CalibPred';
f1val = E_yP_actual;

f2 = 'PredSD';
f2val = V_yP_actual;


mvCalibPred = struct( f1, f1val,...
                    f2, f2val);

end
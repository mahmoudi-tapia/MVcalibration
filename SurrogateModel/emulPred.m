%% Last update: 12/12/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-variate emulator
    

% SAMPLE FUNCTION RECALL 
% [ mvEmulRes ] = emulPred( 'surrogate_model_data.csv', 5, 3, mvEmulRes.rMode , 'surr_pred_data.csv');

function [ mvEmulPred ] = emulPred( emulTrainData, p, q, r_hat, predData)



% Description of the input arguments
% emulTrainData     data filename string - CSV with header
% p                 number of input parameters
% q                 number of output parameters
% r_hat             Vector of size `(p+kappa)x1` which is the estimate of the surrogate model paramaters
% predData          data filename string, Inputs for the purpose of prediction using the emulator without calibration of parameters







% Read the csv files

%'surrogate_model_data.csv'
trainingData = csvread(emulTrainData, 1, 0);  % Ignore the file header

newInput = csvread(predData, 1, 0);  % Ignore the file header

inputRAW = trainingData(:,1:p);    % Collect the input
% Normalize the input -- between 0 and 1
inputNORM = bsxfun(@rdivide,bsxfun(@minus,inputRAW,min(inputRAW)),max(inputRAW)-min(inputRAW));

outputRAW = trainingData(:,p+1:p+q);    % Collect the output matrix
% Normalize output -- mean 0
D = bsxfun(@rdivide,bsxfun(@minus,outputRAW,mean(outputRAW)),std(outputRAW));    

% Generate the design matrix
% Here we take m = 1
% H = designMat (inputNORM, 1);
H = designMat (inputNORM, p+1);



[~, A_valid_chol] = corrMat2(r_hat, inputNORM); % Compute correlation matrix
[ BetaHAT_valid ] = BetaHAT2( A_valid_chol, H, D);
[ SigmaHAT_GLS ] = SigmaHAT2( A_valid_chol, H, D, BetaHAT_valid);


X_RAW = newInput(:,1:p);
% Normalized the input
X_NORM = bsxfun(@rdivide,bsxfun(@minus,X_RAW,min(inputRAW)),max(inputRAW)-min(inputRAW));


% Compute the mean vector (predictions) for the T-distribution
[ mu_valid, T_valid ] = muDstar2( X_NORM, BetaHAT_valid, D, A_valid_chol, H, inputNORM, r_hat);

% Compute the prediction variance
predVAR = cDstar2( X_NORM, T_valid, A_valid_chol, H, r_hat);
% predSD_NORM = sqrt( diag( SigmaHAT_GLS * cDstar( x_validNORM, x_validNORM, A_valid, H, inputNORM, r_hat )));
predSD_NORM = sqrt(bsxfun(@times,diag(predVAR),diag(SigmaHAT_GLS)'));
predSD = bsxfun(@times,predSD_NORM,std(outputRAW));

% Scale back the predictions and compute error
PredictedOutput = bsxfun(@plus,bsxfun(@times,mu_valid,std(outputRAW)),mean(outputRAW));


%% Save the resutls in a struct
f1 = 'EmulPred';
f1val = PredictedOutput;

f2 = 'PredSD';
f2val = predSD;


mvEmulPred = struct( f1, f1val,...
                    f2, f2val);






end
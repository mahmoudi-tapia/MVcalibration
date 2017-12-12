%% Last update: 11/29/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% This script generates the emulator
    

% SAMPLE FUNCTION RECALL 
% [ mvCalibRes ] = mvCalibrator( r_hat, 'surrogate_model_data.csv', 'calibration_data.csv', 5, 3, 2, 500, 4);

function [ mvCalibRes ] = mvCalibrator( r_hat, emulTrainData, calibTrainData, P, q, t, MCMC, k)
tStart = tic;

% Description of the input arguments

% r_hat             Vector of size (p+kappa)x1 which is the estimate of the surrogate model paramaters
% emulTrainData     numerical experiments (simulation) data filename string - CSV with header
% calibTrainData     physical experiments data filename string - CSV with header
% p                 number of input parameters (control variables + calibration parameters)
% q                 number of output parameters
% t                 number of control inputs
% MCMC              number of iterations for MCMC     
% k                 parameter for the leave-k-out cross validation
%                   (Number of data points set aside for validation)  



% ---------------------    
%       PARAMETERS
% ---------------------
% Acquire the data
simData = csvread(emulTrainData, 1, 0);


% sanity check
if ( size(simData,2) ~= P+q)
    disp('The dataset seems to be incorrect. Number of columns in the training data is not equal with the sum of inputs and outputs.')
    return
elseif t >= P
    disp('Number of calibration parameters must be smaller than the total number of input variables.')
    return    
end

% ---------------------
iter = MCMC;   % Number of iterations for MCMC
PARAM2 = 'Truncated Normal - Normal proposal for THETA';    % Note on the proposal distributions
PARAM3 = 'Sim. data size: 130, Calib. data 1 with Tp';  % Note on the experimental data
COMMENT1 = 'Experimental data size 24, with 4 validation chose RANDOMLY ';
COMMENT2 = 'Predictions based on POSTERIOR MEAN';

% Acquire the experimental and prediction data
fullCalibData = csvread(calibTrainData, 1, 0);  % Ignore the file header
Nfull = size(fullCalibData,1);


% sanity check
if (mod(length(fullCalibData),k) ~= 0)
    disp('Bad choice of "k" for leave-k-out cross validation.')
    return
end

% ---------------------
Date=datestr(now,31);
startTime=tic;

% Declare the prediction vector
calibPred = zeros(length(fullCalibData),q);
calibPredSD = zeros(length(fullCalibData),q);

noRun = length(fullCalibData)/k;
% Regular k-fold cross validation
for iRun=1:noRun                
        
        validIndex = [(iRun-1)*k+1:iRun*k]';
        % Run the emulator for the reduced data
        [calibPred(validIndex, :), calibPredSD(validIndex, :)] = ...
        calibratorForCV( iRun, validIndex, Date, startTime, r_hat, simData, fullCalibData, P, q, t, MCMC );
        
end
	

%% Plot the results of cross-validation
cvCalibPlot( q, [fullCalibData(:,t+1:end), calibPred,calibPredSD] );


%% Calibration with the full training data
[ thetaSample, thetaMODE, rDeltaST, psiDeltaST, psiEpsST ] = calibratorAfterCV( Date, startTime, r_hat,...
simData, fullCalibData, P, q, t, MCMC );

%% Save the resutls in a struct
f1 = 'calibPred';
f1val = calibPred;

f2 = 'thetaSample';
f2val = thetaSample;

f3 = 'thetaMode';
f3val = thetaMODE;

f4 = 'rDeltaHAT';
f4val = rDeltaST;

f5 = 'psiDeltaHAT';
f5val = psiDeltaST;

f6 = 'psiEpsHAT';
f6val = psiEpsST;

mvCalibRes = struct( f1, f1val,...
                    f2, f2val,...
                    f3, f3val,...
                    f4, f4val,...
                    f5, f5val,...
                    f6, f6val);

toc(tStart);
end
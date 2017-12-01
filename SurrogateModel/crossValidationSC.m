%% Last update: 09/13/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% This script performs leave-k-out cross validation for an emulator




% ---------------------    
%       TRAINING DATA
% ---------------------

% trainingData = csvread('Data1_rev1.csv', 1, 0);  % Ignore the file header
trainingData = csvread('Data2_extended.csv', 1, 0);  % Ignore the file header
% Data = csvread('Data4_MVGP_alldata_4output.csv', 1, 0);  trainingData = Data(1:60,:);
% trainingData = testSet(60, 'Default');% Use test data set

% ---------------------    
%       PARAMETERS
% ---------------------
p = 5;  % Number of input parameters
q = size(trainingData,2)-p;  % Number of output parameters

PARAM1 = 500;   % Number of iterations for MCMC
PARAM2 = 'Truncated Normal';    % Note on the proposal distribution
PARAM3 = 'FEM data - set 1';  % Note on the training data

% Number of data points set aside for validation
k = 13;
noRun = length(trainingData)/k;

if (mod(length(trainingData),k) ~= 0) 
    disp('Bad choice of "k" for leave-k-out cross validation.')
    return 
end
    
% Declare the prediction vector
prediction = zeros(length(trainingData),q);

% for iRun=1:noRun        
    
    testIndex = [(iRun-1)*k+1:iRun*k]';
    
    % Run the emulator for the reduced data
    prediction(testIndex, :) = emulatorForCrossVal2( iRun, trainingData, testIndex,...
        p, q, PARAM1, PARAM2, PARAM3);
    
% end
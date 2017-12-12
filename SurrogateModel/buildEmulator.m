%% Last update: 11/29/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-variate emulator
    

% SAMPLE FUNCTION RECALL 
% [ mvEmulRes ] = buildEmulator( 'surrogate_model_data.csv', 5, 3, 500, 13);

function [ mvEmulRes ] = buildEmulator( emulTrainData, p, q, MCMC, k)
tStart=tic;

% Description of the input arguments
% emulTrainData     data filename string - CSV with header
% p                 number of input parameters
% q                 number of output parameters
% MCMC              number of iterations for MCMC     
% k                 parameter for k-fold cross validation
%                   (Number of data points set aside for validation)  


%'surrogate_model_data.csv'
trainingData = csvread(emulTrainData, 1, 0);  % Ignore the file header

% sanity check
if ( size(trainingData,2) ~= p+q)
    disp(['The dataset seems to be incorrect. Number of columns in the training data is not equal with the sum of inputs and outputs.'])
    return
end

PARAM1 = MCMC;
PARAM2 = 'Truncated Normal';    % Note on the proposal distribution
PARAM3 = 'FEM data - set 1';  % Note on the training data

% Regular k-fold cross validation
if k>1
    noRun = length(trainingData)/k;
    
    if (mod(length(trainingData),k) ~= 0)
        disp('Bad choice of "k" for leave-k-out cross validation.')
        return
    end
    
    % Declare the prediction vector
    emulPred = zeros(length(trainingData),q);
    emulPredSD = zeros(length(trainingData),q);
    
    for iRun=1:noRun
        
        testIndex = [(iRun-1)*k+1:iRun*k]';
        
        % Run the emulator for the reduced data
        [emulPred(testIndex, :), emulPredSD(testIndex, :)] = ...
            emulatorForCrossVal2( iRun, trainingData, testIndex,...
            p, q, PARAM1, PARAM2, PARAM3);
        
    end
end

%% Plot the results of cross-validation
cvEmulPlot( q, [trainingData(:,p+1:end), emulPred,emulPredSD] );


%% Emulation with the full training data
[ rMode, rSample ] = emulatorAfterCV( trainingData, p, q, PARAM1, PARAM2, PARAM3 );

%% Save the resutls in a struct
f1 = 'EmulPred';
f1val = emulPred;

f2 = 'rSample';
f2val = rSample;

f3 = 'rMode';
f3val = rMode;

mvEmulRes = struct( f1, f1val,...
                    f2, f2val,...
                    f3, f3val);

                

toc(tStart)
end
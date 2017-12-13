# Multivariate statistical calibration of computer models

Please read the following procedure before using the MATLAB codes in this repository.

## Introduction
This repository consists of codes and data files for the purpose of emulation (surrogate modeling) and statistical calibration of computer models using Gaussian processes. The data used here as an example come from an FEM thermal model developed for metal-based additive manufacturing. The emulation and calibration are performed for a multivariate model with the following properties:

- Number of control inputs: `p`
- Number of calibration parameters: `t`
- Number of outputs: `q`

The user can use the repository for building and validation a surrogate model from her own data. This task is explained in **Step 1**.
In addition, the user can provide experimental data to calibrate her model parameters. This task is explained in **Step 2**.

### The data

We use a number of CSV spreadsheets throughout this project, as listed in the following table. Note that the the spreadsheets are assumed to have header in all cases.

Description of the spreadsheet | Sample spreadsheet used in the repository
------------ | -------------
Numerical experiments data, i.e. the inputs, parameters, and outputs from the computer simulation model | `surrogate_model_data.csv`
Inputs for the purpose of prediction using the emulator without calibration of parameters | `surr_pred_data.csv`
Physical experiments data, i.e. the inputs (without parameters) and outputs from experimentation | `calibration_data.csv`
Inputs for the purpose of prediction using the emulator after calibration is done| `calib_pred_data.csv`






## Step 1: Emulation (Surrogate modeling)

The first step is to build an emulator (a surrogate model) using the numerical experiments data. Use the following procedure:

1. Go to the folder "Surrogate Model"
2. Use the function `buildEmulator` with the following inputs:

* **emulTrainData**: 	Data filename string for a spreadsheet with `N` rows and `(p+t+q)` columns
* **p+t**: 		Summation of number of control inputs and calibration parameters
* **q**:		number of output parameters
* **MCMC**: 	number of iterations for MCMC
* **k**: 		parameter for k-fold cross validation, (Number of data points set aside for validation) 


A screenshot of the example numerical data used is shown below:

![screenshot of the data csv](https://github.com/mahmoudi-tapia/MVcalibration/blob/master/SurrogateModel/surrogate_model_data_scr.JPG)

A sample function recall for a model with 2 control inputs, 3 calibration parameters, and 3 outputs will look like:

```javascript
[ mvEmulRes ] = buildEmulator ('surrogate_model_data.csv', 5, 3, 5000, 13);
```
The output **mvEmulRes** will be a struct with the following elements:

* **EmulPred**: 	Predictions for each of the training data input using the surrogate model
* **rSample**: 		Posterior distribution for the roughness parameters for the surrogate model
* **rMode**:		Mode of the posterior distribution for the roughness parameters for the surrogate model

### Prediction without calibrated parameters
Once the emulator is ready, it can be used predictions either without or with the calibrated parameters. A specific function `emulPred` is written for when the user wants to predict using the emulator without the calibrated parameters. In this case, user should specify her desired control inputs in addition to the parameters in a separate spreadsheet. Use the function `emulPred` with the following inputs:

* **emulTrainData**: 	Data filename string for a spreadsheet with `N` rows and `(p+t+q)` columns
* **p+t**: 		Summation of number of control inputs and calibration parameters
* **q**:		number of output parameters
* **r_hat**:     Vector of size `(p+t)x1` which is the estimate of the surrogate model parameters
* **predData**: Data filename string of the inputs for the purpose of prediction using the emulator without calibration of parameters (a spreadsheet with `S` rows and `(p+t)` columns)

A sample function recall for a model with 2 control inputs, 3 calibration parameters, and 3 outputs will look like:

```javascript
[ mvEmulRes ] = emulPred ('surrogate_model_data.csv', 5, 3, mvEmulRes.rMode , 'surr_pred_data.csv');
```

The output **emulPred** will be a struct with the following elements:

* **EmulPred**: 	Predictions for the input using the surrogate model
* **PredSD**: 		Standard error of the predictions

The output of the emulation step can also be used for calibration of parameters as described in **Step 2**.

## Step 2: Calibration

The second step is to calibrate the emulator that was built in the previous step. Note that the output of the previous step is need; particularly, a vector of the estimate of the roughness parameters is needed. By default, we use the posterior mode as the estimates `r_hat`. To accomplish the multivariate calibration, use the following procedure:

1. Go to the folder "Calibration"
2. Use the function `mvCalibrator` with the following inputs:

* **r_hat**:     Vector of size `(p+t)x1` which is the estimate of the surrogate model parameters
*  **emulTrainData**: 	Data filename string for a spreadsheet with `N` rows and `(p+t+q)` columns
*  **calibTrainData**:   Data filename string for a spreadsheet with `M` rows and `(p+q)` columns
* **p+t**: 		Summation of number of control inputs and calibration parameters
* **q**:		number of output parameters
* **t**:             number of control inputs
* **MCMC**: 	number of iterations for MCMC
* **k**: 		parameter for k-fold cross validation, (Number of data points set aside for validation)

A sample function recall will look like:

```javascript
[ mvCalibRes ] = mvCalibrator (r_hat, 'surrogate_model_data.csv', 'calibration_data.csv', 5, 3, 2, 5000, 4);
```
The output **mvCalibRes** will be a struct with the following elements:

* **calibPred**: 	Predictions for each of the predictions done after calibration
* **thetaSample**: 		Posterior distribution for the calibration parameters
* **thetaMode**:		Mode of the posterior distribution for the calibration parameters 
* **rDeltaHAT**: 		Estimate of the roughness parameters of the model outputs
* **psiDeltaHAT**:		Estimate of the variance of the model outputs
* **psiEpsHAT**:		Estimate of the noise variance

### Prediction after parameter calibration
Once the parameters are calibrated and estimates for the model discrepancy and noise parameters are determined, more predictions can be done at new input values. A specific function `calibPred` is written for when the user wants to predict using the emulator with the calibrated parameters. In this case, user should specify her desired control inputs in a separate spreadsheet. Use the function `calibPred` with the following inputs:

* **r_hat**:     Vector of size `(p+t)x1` which is the estimate of the surrogate model parameters
*  **emulTrainData**: 	Data filename string for a spreadsheet with `N` rows and `(p+t+q)` columns
*  **calibTrainData**:   Data filename string for a spreadsheet with `M` rows and `(p+q)` columns
* **p+t**: 		Summation of number of control inputs and calibration parameters
* **q**:		number of output parameters
* **t**:             number of control inputs
* **thetaMode**:		Mode of the posterior distribution for the calibration parameters 
* **rDeltaHAT**: 		Estimate of the roughness parameters of the model outputs
* **psiDeltaHAT**:		Estimate of the variance of the model outputs
* **psiEpsHAT**:		Estimate of the noise variance
* **predDataCSV**:       data filename string, Inputs for the purpose of prediction using the emulator after calibrating parameters (a spreadsheet with `S` rows and `p` columns)

A sample function recall for a model with 2 control inputs, 3 calibration parameters, and 3 outputs will look like:

```javascript
[ mvCalibPred ] = calibPred( r_hat, 'surrogate_model_data.csv', 'calibration_data.csv', 5, 3, 2, mvCalibRes.thetaMode , mvCalibRes.rDeltaHAT, mvCalibRes.psiDeltaHAT, mvCalibRes.psiEpsHAT, 'calib_pred_data.csv');
```

The output **mvCalibPred** will be a struct with the following elements:

* **CalibPred**: 	Predictions for the input using the calibrated surrogate model
* **PredSD**: 		Standard error of the predictions

## Remarks

1. For each step (Emulation and Calibration), results of the each fold of the leave-k-out validation is saved as plots and MATLAB workspace files. These files are stored in the folders names 'CodeOutput'.

2. After the cross-validation for each step, another final training is performed with the full data i.e. no data points left out for any validation. The workspace of this final run is also saved in the folder 'CodeOutput' for future reference.

3. In case of issues/questions, please do not hesitate to contact the developers at mahmoudi@tamu.edu. 

 
Good luck! :+1:

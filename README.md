# Multivariate statistical calibration of computer models

Please read the following procedure before using the MATLAB codes in this repository.

### The data

There are two spreadsheets used throughout for this project:

1. surrogate_model_data.csv
2. calibration_data.csv

(https://github.com/mahmoudi-tapia/MVcalibration/blob/master/SurrogateModel/surr_data_scr.JPG)

The first spreadsheet containts the numerical experiments data, i.e. the inputs and outputs from the simulation model. The second spreadsheet contains the physical experiments data. The CSV format with header is used for both.

### Step 1: Emulation (Surrogate modeling)

The first step is to build an emulator (a surrogate model) using the numerical experiments data. Use the following procedure:

1. Go to the folder "Surrogate Model"
2. Use the function `mvEmulator` with the following inputs:

* **emulTrainData**: 	Data filename string - CSV with header
* **p**: 		number of input parameters
* **q**:		number of output parameters
* **MCMC**: 	number of iterations for MCMC
* **k**: 		parameter for k-fold cross validation, (Number of data points set aside for validation) 

A sample function recall will look like:

```javascript
[ mvEmulRes ] = mvEmulator ( 'surrogate_model_data.csv', 5, 3, 500, 13);
```
The output **mvEmulRes** will be a struct with the following elements:

* **EmulPred**: 	Predictions for each of the training data input using the surrogate model
* **rSample**: 		Posterior distribution for the roughness parameters for the surrogate model
* **rMode**:		Mode of the posterior distribution for the roughness parameters for the surrogate model

The output of the emulation step can then be used for the calibration described in Step 2.

### Step 2: Calibration

The second step is to calibrate the emulator that was built in the previous step. Note that the output of the previous step is need; particularly, a vector of the estiamte of the roughness parameters is needed. By default, we use the posterior mode as the estimates `r_hat`. To accomplish the multivariate calibration, use the following procedure:

1. Go to the folder "Calibration"
2. Use the function `mvCalibrator` with the following outputs:

* **r_hat**:     Vector of size `(p+kappa)x1` which is the estimate of the surrogate model paramaters
*  **emulTrainData**: 	Data filename string - CSV with header
*  **calibTrainData**:     physical experiments data filename string - CSV with header
* **p**: 		number of input parameters
* **q**:		number of output parameters
* **kappa**:             number of control variables
* **MCMC**: 	number of iterations for MCMC
* **k**: 		parameter for k-fold cross validation, (Number of data points set aside for validation)

A sample function recall will look like:

```javascript
[ mvCalibRes ] = mvCalibrator( r_hat, 'surrogate_model_data.csv', 'calibration_data.csv', 5, 3, 2, 500, 4);
```
The output **mvCalibRes** will be a struct with the following elements:

* **calibPred**: 	Predictions for each of the predictions done after calibration
* **thetaSample**: 		Posterior distribution for the calibration parameters
* **thetaMode**:		Mode of the posterior distribution for the calibration parameters 

### Remarks

1. For each step (Emulation and Calibration), results of the each fold of the leave-k-out validation is saved as plots and MATLAB workspace files. These files are stored in the folders names 'CodeOutput'.

2. After the cross-validation for each step, another final training is performed with the full data i.e. no data points left out for any validation. The workspace of this final run is also saved in the folder 'CodeOutput' for future referece.

3. In case of issues/questions, please do not hesitate to contact the developers at mahmoudi@tamu.edu. 

 
Good luck! :+1:

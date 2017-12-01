%% Last update: 06/11/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Calibrator

% This function gives the logarithm of the prior \pi_{\boldsymbol{R}}(\boldsymbol{r}_delta) for
% the roughness parameters r_\delta
% It is a normalized scale of the product of i.i.d. vague (albeit proper) Log-Logistic prior type 
% LaTeX formula: $\pi_{\boldsymbol{R}}(\boldsymbol{r})=\prod_{i=1}^p (1+r_i^2)^{-1}$

function [ logP ] = pi_r_delta( r )

temp = r.^2;
% logP = logP+1;
logP = -sum(log(1+temp));

end
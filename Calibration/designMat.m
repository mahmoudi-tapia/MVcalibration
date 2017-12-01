%% Last update: 02/13/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% This function gives the design matrix
% This is an incomplete code and only considers m=1 or m=p+1

function [ H ] = designMat( input, m )

% Case m = 1
if m==1
    H = ones(size(input,1),1);
end

% Case m = p+1
if m == size(input,2)+1
    H = [ones(size(input,1),1), input];
end
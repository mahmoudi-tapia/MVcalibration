function [ z ] = propExpPDF( X, mu )
%% Last update: 02/18/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% Multi-variate exponential probability density function

% Sanity check
if length(X) ~= length(mu)
    disp('BAD INPUT! -- CHECK THE CODE')
    return
end

p = length(X);
prob = zeros(p,1);

for i=1:p
    pd = makedist('Exponential', 'mu', mu(i));
    prob(i)=pdf(pd, X(i));
end

z =  prod(prob);


end

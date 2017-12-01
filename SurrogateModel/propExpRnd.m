function [ z ] = propExpRnd( mu )
%% Last update: 02/18/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% Multi-variate exponential random number generator

p = length(mu);
z = zeros(p,1);

for i=1:p
    pd = makedist('Exponential', 'mu', mu(i));
    z(i)=random(pd);
end



end


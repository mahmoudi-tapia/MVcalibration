function [ z ] = propNormRnd( mu, sig, L, H)
%% Last update: 02/22/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% Truncated multi-variate independent normal random number generator

p = length(mu);
z = zeros(p,1);

for i=1:p
    pd = makedist('Normal', 'mu', mu(i),'sigma',sig);
    t = truncate(pd,L,H);
    z(i)=random(t);
end



end
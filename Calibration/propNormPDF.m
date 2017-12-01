function [ z ] = propNormPDF( X, mu, sig, L, H )
%% Last update: 02/18/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% Truncated multi-variate independent normal probability density function

% Sanity check
if length(X) ~= length(mu)
    disp('BAD INPUT! -- CHECK THE CODE')
    return
end

p = length(X);
prob = zeros(p,1);

for i=1:p   
    pd = makedist('Normal', 'mu', mu(i),'sigma',sig);
    t = truncate(pd,L,H);
    prob(i)=pdf(t, X(i));
end

z =  prod(prob);


end

% pd = makedist('Normal',u,v);
% t = truncate(pd,l,h);
% p = pdf(t,x);
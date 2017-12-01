function r = tnormrnd(u,v,l,h,varargin)
% Truncated normal distribution random number generator

pd = makedist('Normal',u,v);
t = truncate(pd,l,h);
if nargin > 4
    m = varargin{1};
    n = varargin{2};
else
    m = 1;
    n = 1;
end
r = random(t,m,n);
return
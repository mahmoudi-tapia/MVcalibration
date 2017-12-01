function p = tnormpdf(x,u,v,l,h)
% Truncated normal distribution probability distribution function

pd = makedist('Normal',u,v);
t = truncate(pd,l,h);
p = pdf(t,x);
return
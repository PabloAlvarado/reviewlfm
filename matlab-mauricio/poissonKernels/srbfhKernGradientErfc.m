function g = srbfhKernGradientErfc(sigmax, lengthX, s1, s2, w, gamma, n, covGrad, erfcz1, erfcz2)

% SRBFHKERNGRADIENT Gradient of the parameters of a SRBFH kernel.
% FORMAT
% DESC computes the gradient of a SRBFH kernel parameter.
% ARG sigmax : length-scale of the spatial gp prior.
% ARG lengthX : length of the spatial domain
% ARG s1 : row inputs for which kernel is to be computed. 
% ARG s2 : column inputs for which kernel is to be computed. 
% ARG w  : precomputed constant.
% ARG gamma : precomputed constant.
% ARG n : integer indicating first series term
% ARG covGrad : partial derivatives
% ARG erfcz1 : precomputed factors
% ARG erfcz2 : precomputed factors
% RETURN g : gradient of the parameters.
%
% SEEALSO : multiKernParamInit, multiKernCompute, srbfhKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2016

% KERN

sinS1 = sin(w(n)*s1);
bTerm = sigmax*gamma(n)/2;
argz2 = (s2-lengthX)/sigmax;
argz1 =  s2/sigmax;
z1 = argz1 + bTerm;
z2 = argz2 + bTerm;

if nargin < 9    
    erfcz1 = Faddeeva_erfc(z1);
    erfcz2 = Faddeeva_erfc(z2);    
end

W = exp(bTerm.^2)*exp(gamma(n)*s2)*(erfcz2 - erfcz1);

dWdsx = ((gamma(n)^2)*sigmax/2)*W - ...
    (2/sqrt(pi))*exp(bTerm.^2)*exp(gamma(n)*s2)*...
    (exp(-(z2^2))*(-((s2-lengthX)/(sigmax^2)) + gamma(n)/2) - ...
     exp(-(z1^2))*(-s2/(sigmax^2) + gamma(n)/2));

dvec = (sqrt(pi)/2)*imag(W) + (sigmax*sqrt(pi)/2)*imag(dWdsx);

g = sum(sum((sinS1*dvec').*covGrad));



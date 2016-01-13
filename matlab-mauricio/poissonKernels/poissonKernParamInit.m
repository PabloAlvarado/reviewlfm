function kern = poissonKernParamInit(kern)

% POISSONKERNPARAMINIT POISSON kernel parameter initialisation.
% The Poisson kernel corresponds to the covariance function for the output
% process of a partial differential equation that corresponds to the
% Poisson equation in two spatial dimensions
%
% \frac{\partial^2 v(x,y)}{\partial x^2} + 
% \frac{\partial^2 v(x,y)}{\partial y^2} = Sf(t) 
%
% where S is the sensitivity coefficient. The kernel also contains the 
% inverse widths associated to RBF covariance of the Gaussian process prior 
% imposed over f(x,y), one for each spatial variable. By default, the
% sensitivity is negative. The particular kernel derived from this equation 
% assumes that the spatial interval is fixed and it can be given as an 
% option. Otherwise, it is assumed to be 1 by default.
%
% Also, the solution to this partial differential equation is given in the
% form of a series. We also specify the number of terms in the series. By
% default the number of terms is 5.
%
% FORMAT
% DESC initialises the heat kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2016

% KERN

if kern.inputDimension ~= 2
  error('POISSON kernel only valid for 2-D inputs.')
end

if isfield(kern, 'options') && isfield(kern.options, 'nTerms')
    kern.nTerms = kern.options.nTerms;
else
    kern.nTerms = 20;
end

if isfield(kern, 'options') && isfield(kern.options, 'lengthX')
    kern.lengthX = kern.options.lengthX;
else
    kern.lengthX = 1;
end

if isfield(kern, 'options') && isfield(kern.options, 'lengthY')
    kern.lengthY = kern.options.lengthY;
else
    kern.lengthY = 1;
end

kern.inverseWidthX = 1;
kern.inverseWidthY = 1;
kern.sensitivity = 1;
kern.nParams = 3;
kern.transforms.type = optimiDefaultConstraint('positive');
kern.transforms.index = 1:(kern.nParams-1); % The sensitivity could be negative
kern.isStationary = false; 






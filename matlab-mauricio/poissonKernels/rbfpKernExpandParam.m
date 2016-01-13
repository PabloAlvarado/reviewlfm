function kern = rbfpKernExpandParam(kern, params)

% RBFPKERNEXPANDPARAM Create kernel structure from RBFP kernel's parameters.
% FORMAT
% DESC returns a radial basis function POISSON kernel structure filled with 
% the parameters in the given vector. This is used as a helper function to
% enable parameters to be optimised in, for example, the NETLAB
% optimisation functions.
% ARG kern : the kernel structure in which the parameters are to be
% placed.
% ARG param : vector of parameters which are to be placed in the
% kernel structure.
% RETURN kern : kernel structure with the given parameters in the
% relevant locations.
%
% SEEALSO : rbfpKernParamInit, rbfpKernExtractParam, kernExpandParam
%
% COPYRIGHT : Mauricio A. Alvarez, 2016

% KERN

kern.inverseWidthX = params(1);
kern.inverseWidthY = params(2);

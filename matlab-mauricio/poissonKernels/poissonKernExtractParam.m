function [params, names] = poissonKernExtractParam(kern)

% POISSONKERNEXTRACTPARAM Extract parameters from the Poisson kernel structure.
% FORMAT
% DESC Extract parameters from the Poisson kernel structure into a vector of
% parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC Extract parameters and their names from the Poisson kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO poissonKernParamInit, poissonKernExpandParam, scg, conjgrad
%
% COPYRIGHT : Mauricio A. Alvarez, 2016
%
% KERN



params = [kern.inverseWidthX kern.inverseWidthY kern.sensitivity];
if nargout > 1
    names = {'inverse width space X.', 'inverse width space Y.', 'sensitivity'};
end

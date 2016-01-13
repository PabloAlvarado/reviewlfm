function k = rbfpKernDiagCompute(kern, x)

% RBFPKERNDIAGCOMPUTE Compute diagonal of RBFP kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the radial basis 
% function Poisson kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : rbfpKernParamInit, kernDiagCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2016

% KERN

% Since the variance is on the diagonal is just ones

if size(x, 2) ~= 2
    error('Input can only have two columns');
end

k = ones(size(x,1),1);

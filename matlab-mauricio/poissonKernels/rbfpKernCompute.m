function K = rbfpKernCompute(kern, x1, x2)

% RBFPKERNCOMPUTE Compute the RBFP kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the radial basis function Poisson
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x1 : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN K : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the radial basis function Poisson
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x1 : input data matrix in the form of a design matrix.
% RETURN K : the kernel matrix computed at the given points.
%
% SEEALSO : rbfpKernParamInit, kernCompute, 
%
% COPYRIGHT : Mauricio A. Alvarez, 2016

% KERN

if nargin < 3
    x2 = x1;
end

if size(x1, 2) ~= 2 || size(x2, 2) ~= 2
    error('Input can only have two columns');
end


% Split the domain into the x spatial domain and y spatial domain

sx1 = x1(:,1);
sx2 = x2(:,1);
sy1 = x1(:,2);
sy2 = x2(:,2);

kern.rbf.inverseWidth = kern.inverseWidthX;
Kx = rbfKernCompute(kern.rbf, sx1, sx2);
kern.rbf.inverseWidth = kern.inverseWidthY;
Ky = rbfKernCompute(kern.rbf, sy1, sy2);
 
K = Kx.*Ky;




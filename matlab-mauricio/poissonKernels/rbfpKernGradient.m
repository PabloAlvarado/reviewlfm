function g = rbfpKernGradient(kern, x1, x2, covGrad)

% RBFPKERNGRADIENT Gradient of RBFP kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the
% radial basis function Poisson kernel's parameters. As well as the kernel
% structure and the input positions, the user provides a matrix PARTIAL
% which gives the partial derivatives of the function with respect to the
% relevant elements of the kernel matrix.
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x1 : the input locations for which the gradients are being
% computed.
% ARG covGrad : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The argument takes
% the form of a square matrix of dimension  numData, where numData is
% the number of rows in X.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters. The ordering of the vector should match
% that provided by the function kernExtractParam.
%
% FORMAT
% DESC computes the derivatives as above, but input locations are
% now provided in two matrices associated with rows and columns of
% the kernel matrix.
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x1 : the input locations associated with the rows of the
% kernel matrix.
% ARG x2 : the input locations associated with the columns of the
% kernel matrix.
% ARG covGrad : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The matrix should
% have the same number of rows as X1 and the same number of columns
% as X2 has rows.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters.
%
% SEEALSO: rbfpKernParamInit, kernGradient
%
% COPYRIGHT : Mauricio A. Alvarez, 2016

% KERN

if nargin < 4
    covGrad = x2;
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

covGradx = covGrad.*Ky;
covGrady = covGrad.*Kx;

kern.rbf.inverseWidth = kern.inverseWidthX;
gx = rbfKernGradient(kern.rbf, sx1, sx2, covGradx);
kern.rbf.inverseWidth = kern.inverseWidthY;
gy = rbfKernGradient(kern.rbf, sy1, sy2, covGrady);

g = [gx(1) gy(1)];


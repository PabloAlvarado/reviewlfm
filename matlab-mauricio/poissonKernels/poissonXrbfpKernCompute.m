function [K, sK] = poissonXrbfpKernCompute(poissKern, rbfpKern, x1, x2)

% POISSONXRBFPKERNCOMPUTE Cross kernel between a POISSON and a RBFP kernels.
% FORMAT
% DESC computes cross kernel terms between a POISSON kernel and a RBFP kernel
% for the multiple output kernel.
% ARG poissKern : the kernel structure associated with the POISSON kernel.
% ARG rbfpKern : the kernel structure associated with the RBFP kernel.
% ARG x1 : inputs for which kernel is to be computed. First column represent
% the time points, while the second column represents the spatial points.
% Entries with Inf indicate missing values.
% RETURN K : block of values from kernel matrix.
% RETURN sK : unscaled kernel matrix
%
% FORMAT
% DESC computes cross kernel terms between a POISSON kernel and a RBFP kernel
% for the multiple output kernel.
% ARG poissKern : the kernel structure associated with the POISSON kernel.
% ARG rbfpKern : the kernel structure associated with the RBFP kernel.
% ARG x1 : row inputs for which kernel is to be computed. First column
% corresponds to time points and the second column corresponds to spatial
% points. Entries with Inf indicate missing values.
% ARG x2 : column inputs for which kernel is to be computed. First column
% corresponds to time points and the second column corresponds to spatial
% points. Entries with Inf indicate missing values.
% RETURN k : block of values from kernel matrix.
% RETURN sK : unscaled kernel matrix
%
% SEEALSO : multiKernParamInit, multiKernCompute, poissonKernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2016

% KERN

if nargin < 4
    x2 = x1;
end
if size(x1, 2) ~= 2 || size(x2, 2) ~= 2
    error('Input can only have two columns');
end
if (poissKern.inverseWidthX ~= rbfpKern.inverseWidthX) || ...
        (poissKern.inverseWidthY ~= rbfpKern.inverseWidthY)
    error('Kernels cannot be cross combined if they have different inverse widths.')
end

% Split the domain into the x spatial domain and y spatial domain

sx1 = x1(:,1);
sx2 = x2(:,1);
sy1 = x1(:,2);
sy2 = x2(:,2);

K = zeros(length(sx1), length(sx2));

% Although this is done in heatKernExpandParam.m, we do it here again as a
% precaution.

rbfpKern.rbf.inverseWidth = rbfpKern.inverseWidthX;

sigmax = sqrt(2/poissKern.inverseWidthX);
sigmay = sqrt(2/poissKern.inverseWidthY);
lengthX = poissKern.lengthX;
lengthY = poissKern.lengthY;
nterms = poissKern.nTerms;

% Precompute some terms for spatial variable X
pn = ((1:nterms)*(pi/lengthX))';
gammapn = sqrt(-1)*pn;

% Precompute some terms for spatial variable Y
qm = ((1:nterms)*(pi/lengthY))';
gammaqm = sqrt(-1)*qm;

% Precompute some terms
cK = 4/(lengthX*lengthY);

for n=1:nterms
    Kx = srbfhKernCompute(sigmax, lengthX, sx1, sx2, pn, gammapn, n);
    for m=1:nterms        
        Ky = srbfhKernCompute(sigmay, lengthY, sy1, sy2, qm, gammaqm, m);
        pn2qm2 = pn(n)^2 + qm(m)^2;
        K = K + (Kx.*Ky)/(pn2qm2);
    end
end
sK = cK*K;
K = poissKern.sensitivity*sK;





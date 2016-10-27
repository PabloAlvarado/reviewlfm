function [g1, g2] = poissonXrbfpKernGradient(poissKern, rbfpKern, x1, x2, covGrad)

% POISSONXRBFPKERNGRADIENT Gradient wrt parameters between a POISSON and a RBFP.
% FORMAT
% DESC computes the gradients wrt parameters of a cross kernel term between
% a POISSON kernel and a RBFP kernel for the multiple output kernel.
% ARG poissKern : the kernel structure associated with the POISSON kernel.
% ARG rbfpKern : the kernel structure associated with the RBFP kernel.
% ARG x1 : inputs for which kernel is to be computed. First column represent
% the time points, while the second column represents the spatial points.
% Entries with Inf indicate missing values.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see poissonKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see rbfKernExtractParam.
%
% FORMAT
% DESC computes the gradients wrt parameters of a cross kernel term between
% a POISSON kernel and a RBFP kernel for the multiple output kernel.
% ARG poissKern : the kernel structure associated with the POISSON kernel.
% ARG rbfpKern : the kernel structure associated with the RBFP kernel.
% ARG x1 : row inputs for which kernel is to be computed. First column
% corresponds to time points and the second column corresponds to spatial
% points. Entries with Inf indicate missing values.
% ARG x2 : column inputs for which kernel is to be computed. First column
% corresponds to time points and the second column corresponds to spatial
% points. Entries with Inf indicate missing values.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see poissonKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see rbfpKernExtractParam.
%
% SEEALSO : poissonXrbfpKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2016

% KERN

if nargin < 5
    covGrad = x2;
    x2 = x1;
end
if size(x1, 2) ~= 2 || size(x2, 2) ~= 2
    error('Input can only have two columns');
end
if (poissKern.inverseWidthX ~= rbfpKern.inverseWidthX) || ...
        (poissKern.inverseWidthY ~= rbfpKern.inverseWidthY)
    error('Kernels cannot be cross combined if they have different inverse widths.')
end

sx1 = x1(:,1);
sx2 = x2(:,1);
sy1 = x1(:,2);
sy2 = x2(:,2);

sK = zeros(length(sx1), length(sx2));

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

g1 = zeros(1,3);
g2 = zeros(1,2);

%%%%% Code with Faddeeva_w
% for n=1:nterms
%     [Kx, wzx1, wzx2] = srbfhKernCompute(sigmax, lengthX, sx1, sx2, ...
%             pn, gammapn, n);
%     covGrady = covGrad.*Kx;        
%     for m=1:nterms
%         [Ky, wzy1, wzy2] = srbfhKernCompute(sigmay, lengthY, sy1, sy2, ...
%             qm, gammaqm, m);        
%         covGradx = covGrad.*Ky;        
%         pn2qm2 = pn(n)^2 + qm(m)^2;
%         gx = srbfhKernGradient(sigmax, lengthX, sx1, sx2, pn, ...
%             gammapn, n, covGradx, wzx1, wzx2);        
%         gx = gx/pn2qm2;
%         gx = -(1/sqrt(2*poissKern.inverseWidthX^3))*gx; % Transforms to the derivative of the inverse width
%         g1(1) = g1(1) + gx;
%         gy = srbfhKernGradient(sigmay, lengthY, sy1, sy2, qm, ...
%             gammaqm, m, covGrady, wzy1, wzy2);
%         gy = gy/pn2qm2;
%         gy = -(1/sqrt(2*poissKern.inverseWidthY^3))*gy; % Transforms to the derivative of the inverse width
%         g1(2) = g1(2) + gy;
%         sK = sK + (Kx.*Ky)/pn2qm2;
%     end
% end

%%%%% Code with Faddeeva_erfc
for n=1:nterms
    [Kx, erfcx1, erfcx2] = srbfhKernComputeErfc(sigmax, lengthX, sx1, sx2, ...
            pn, gammapn, n);
    covGrady = covGrad.*Kx;        
    for m=1:nterms
        [Ky, erfcy1, erfcy2] = srbfhKernComputeErfc(sigmay, lengthY, sy1, sy2, ...
            qm, gammaqm, m);        
        covGradx = covGrad.*Ky;        
        pn2qm2 = pn(n)^2 + qm(m)^2;
        gx = srbfhKernGradientErfc(sigmax, lengthX, sx1, sx2, pn, ...
            gammapn, n, covGradx, erfcx1, erfcx2);        
        gx = gx/pn2qm2;
        gx = -(1/sqrt(2*poissKern.inverseWidthX^3))*gx; % Transforms to the derivative of the inverse width
        g1(1) = g1(1) + gx;
        gy = srbfhKernGradientErfc(sigmay, lengthY, sy1, sy2, qm, ...
            gammaqm, m, covGrady, erfcy1, erfcy2);
        gy = gy/pn2qm2;
        gy = -(1/sqrt(2*poissKern.inverseWidthY^3))*gy; % Transforms to the derivative of the inverse width
        g1(2) = g1(2) + gy;
        sK = sK + (Kx.*Ky)/pn2qm2;
    end
end

g1(1:2) = poissKern.sensitivity*cK*g1(1:2);
g1(3) = cK*(sum(sum(covGrad.*sK)));






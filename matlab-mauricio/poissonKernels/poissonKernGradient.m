function g = poissonKernGradient(poissKern, x, varargin)

% POISSONKERNGRADIENT Gradient of POISSON kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the Poisson
% kernel's parameters. As well as the kernel structure and the
% input positions, the user provides a matrix PARTIAL which gives
% the partial derivatives of the function with respect to the
% relevant elements of the kernel matrix.
% ARG poissKern : the kernel structure for which the gradients are being
% computed.
% ARG x : the input locations for which the gradients are being
% computed.
% ARG covGrad : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The argument takes
% the form of a square matrix of dimension  numData, where numData is
% the number of rows in X.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters. The ordering of the vector should match
% that provided by the function kernExtractParam.
%
% SEEALSO poissonKernParamInit, kernGradient, poissonKernDiagGradient, kernGradX
%
% COPYRIGHT : Mauricio A. Alvarez, 2016

% KERN

if size(x, 2) ~= 2 
    error('Input can only have two columns');
end

if nargin < 4
    x2 = x;
    covGrad = varargin{1};
else
    x2 = varargin{1};
    covGrad = varargin{2};
end

% Split the domain into the x spatial domain and y spatial domain

sx1 = x(:,1);
sy1 = x(:,2);
sx2 = x2(:,1);
sy2 = x2(:,2);

sK = zeros(length(sx1), length(sx2));

sigmax = sqrt(2/poissKern.inverseWidthX);
sigmay = sqrt(2/poissKern.inverseWidthY);
lengthX = poissKern.lengthX;
lengthY = poissKern.lengthY;
nterms = poissKern.nTerms;

% Precompute some terms for spatial variable X
pn = ((1:nterms)*(pi/lengthX))';
gammapn = sqrt(-1)*pn;
z1gpn = sigmax*gammapn/2;
z2gpn = lengthX/sigmax + z1gpn;
wz1gpn = wofzPoppe(sqrt(-1)*z1gpn);
wz2gpn = wofzPoppe(sqrt(-1)*z2gpn);

% Precompute some terms for spatial variable Y
qm = ((1:nterms)*(pi/lengthY))';
gammaqm = sqrt(-1)*qm;
z1gqm = sigmay*gammaqm/2;
z2gqm = lengthY/sigmay + z1gqm;
wz1gqm = wofzPoppe(sqrt(-1)*z1gqm);
wz2gqm = wofzPoppe(sqrt(-1)*z2gqm);

% Constant tmer in front of the kernel

cK = 16/((lengthX^lengthY)^2);
g = zeros(1, poissKern.nParams);

for n=1:nterms
    for np=1:nterms
        if  (mod(n+np,2)==0)
            Kx = sheatKernCompute(sigmax, lengthX, sx1, sx2, pn, ...
                gammapn, wz1gpn, wz2gpn, n, np);
            for m=1:nterms
                for mp=1:nterms
                    if  (mod(m+mp,2)==0)
                        pn2qm2 = pn(n)^2 + qm(m)^2;
                        pnp2qmp2 = pn(np)^2 + qm(mp)^2;
                        Ky = sheatKernCompute(sigmay, lengthY, sy1, sy2, qm, ...
                            gammaqm, wz1gqm, wz2gqm, m, mp);
                        % Derivative wrt sigmax   
                        covGradx = covGrad.*Ky;
                        gx = sheatKernGradient(sigmax, lengthX, sx1, sx2, pn, ...
                            gammapn, wz1gpn, wz2gpn, n, np, covGradx);
                        gx = gx/(pn2qm2*pnp2qmp2);
                        gx = -(1/sqrt(2*poissKern.inverseWidthX^3))*gx;
                        g(1) = g(1) + gx;
                        % Derivative wrt sigmay                        
                        covGrady = covGrad.*Kx;
                        gy = sheatKernGradient(sigmay, lengthY, sy1, sy2, qm, ...
                            gammaqm, wz1gqm, wz2gqm, m, mp, covGrady);
                        gy = gy/(pn2qm2*pnp2qmp2);
                        gy = -(1/sqrt(2*poissKern.inverseWidthY^3))*gy;
                        g(2) = g(2) + gy;
                        
                        sK = sK + (Kx.*Ky)/(pn2qm2*pnp2qmp2);
                    end
                end
            end
        end
    end
end

g(1:2) = cK*(poissKern.sensitivity^2)*g(1:2);
g(3) = 2*poissKern.sensitivity*cK*sum(sum(sK.*covGrad));

end





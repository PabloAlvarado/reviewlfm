function g = poissonKernGradient(poissKern, x, covGrad)

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


% Split the domain into the x spatial domain and y spatial domain

sx = x(:,1);
sy = x(:,2);

K = zeros(length(sx));

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

cK = 16/(lengthX^lengthY);
g = zeros(1, poissKern.nParams);

for n=1:nterms
    for m=1:nterms
        for np=1:nterms
            for mp=1:nterms
                if (mod(n+np,2)==0) && (mod(m+mp,2)==0)
                    %%%%%
                    epsilon = 1e-9;
                    param = sigmax + epsilon;
                    f1 = computeCvv(param, lengthX, gammapn, wz1gpn, wz2gpn, n, np);
                    param = sigmax - epsilon;
                    f2 = computeCvv(param, lengthX, gammapn, wz1gpn, wz2gpn, n, np);
                    dCvvx = 0.5*(f1 - f2)/epsilon;
                    %param = sigmax;
                    %dCvvx = gradientCvv(sigmax, lengthX, gammapn, wz1gpn, wz2gpn, n, np);
%                     
                    %%%%%
                    Cvvx = computeCvv(sigmax, lengthX, gammapn, wz1gpn, wz2gpn, n, np);
                    Cvvy = computeCvv(sigmay, lengthY, gammaqm, wz1gqm, wz2gqm, m, mp);
                    %dCvvx = gradientCvv(sigmax, lengthX, gammapn, wz1gpn, wz2gpn, n, np);
                    dCvvy = gradientCvv(sigmay, lengthY, gammaqm, wz1gqm, wz2gqm, m, mp);
                    pn2qm2 = pn(n)^2 + qm(m)^2;
                    pnp2qmp2 = pn(np)^2 + qm(mp)^2;
                    const = (Cvvx*Cvvy)/(pn2qm2*pnp2qmp2);
                    sinx = sin(pn(n)*sx);
                    csinxsiny = (const*sinx).*sin(qm(m)*sy);
                    sinxsiny = sinx.*sin(qm(m)*sy);
                    sinxpsinyp = sin(pn(np)*sx).*sin(qm(mp)*sy);
                    K = K + csinxsiny*sinxpsinyp';
                    Kterm = sinxsiny*sinxpsinyp';
                    const1 = sum(sum(Kterm.*covGrad));
                    gC = (dCvvx*Cvvy*const1)/(pn2qm2*pnp2qmp2);
                    gC = -(1/sqrt(2*poissKern.inverseWidthX^3))*gC;
                    g(1) = g(1) + gC;
                    gC = (dCvvy*Cvvx*const1)/(pn2qm2*pnp2qmp2);
                    gC = -(1/sqrt(2*poissKern.inverseWidthX^3))*gC;
                    g(2) = g(2) + gC;
                end
            end
        end
    end
end

g(1:2) = cK*(poissKern.sensitivity^2)*g(1:2);
g(3) = 2*poissKern.sensitivity*cK*sum(sum(K.*covGrad));

end





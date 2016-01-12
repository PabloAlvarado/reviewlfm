function [k, sk] = poissonKernDiagCompute(poissKern, x)

% POISSONKERNDIAGCOMPUTE Diagonal of a kernel matrix for a Poisson kernel.
% DESC computes the diagonal of the kernel matrix for the Poisson kernel
% given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
% RETURN sk : unscaled version of the diagonal.
% RETURN sk : unscaled version of the diagonal for the kernel matrix of the
% initial conditions.
%
% SEEALSO : poissonKernParamInit, kernDiagCompute, kernCreate, poissonKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez, 2016

% KERN

if size(x, 2) ~= 2
    error('Input can only have two columns');
end

% Split the domain into the x spatial domain and y spatial domain

sx1 = x(:,1);
sy1 = x(:,2);

k = zeros(length(sx1),1);

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

ck = 16/(lengthX^lengthY);

for n=1:nterms
    for m=1:nterms
        for np=1:nterms
            for mp=1:nterms
                if (mod(n+np,2)==0) && (mod(m+mp,2)==0)
                    Cvvx = computeCvv(sigmax, lengthX, gammapn, wz1gpn, wz2gpn, n, np);
                    Cvvy = computeCvv(sigmay, lengthY, gammaqm, wz1gqm, wz2gqm, m, mp);
                    pn2qm2 = pn(n)^2 + qm(m)^2;
                    pnp2qmp2 = pn(np)^2 + qm(mp)^2;
                    const = (Cvvx*Cvvy)/(pn2qm2*pnp2qmp2);                    
                    csinxsiny = const*sin(pn(n)*sx1);
                    sinxsiny = csinxsiny.*sin(qm(m)*sy1);
                    sinxpsinyp = sin(pn(np)*sx1).*sin(qm(mp)*sy1);                    
                    k = k + sinxsiny.*sinxpsinyp;
                end
            end
        end
    end
end

sk = ck*k;
sk = real(sk);
k = (poissKern.sensitivity^2)*sk;
k = real(k);


end


function [K, sK] = poissonKernCompute(poissKern, x1, x2)

% HEATKERNCOMPUTE Compute a kernel matrix for a Poisson kernel.
% FORMAT
% DESC computes the kernel matrix for the Poisson kernel function given
% inputs associated with rows and columns.
% ARG poissKern : the kernel structure associated with the Poisson kernel
% ARG x1 : inputs for which kernel is to be computed. First column represent
% the time points, while the second column represents the spatial points.
% Entries with Inf indicate missing values.
% RETURN K : block of values from kernel matrix.
%
% FORMAT
% FORMAT
% DESC computes the kernel matrix for the Poisson kernel function given
% inputs associated with rows and columns.
% ARG poissKern : the kernel structure associated with the Poisson kernel
% ARG x1 : inputs for which kernel is to be computed. First column represent
% the time points, while the second column represents the spatial points.
% Entries with Inf indicate missing values.
% RETURN K : block of values from kernel matrix.
% RETURN sK : unscaled kernel matrix
%
% SEEALSO : multiKernParamInit, multiKernCompute, poissonKernParamInit
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

K = zeros(length(sx1), length(sx2));

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


for n=1:nterms
    for np=1:nterms
        if (mod(n+np,2)==0)
            Kx = sheatKernCompute(sigmax, lengthX, sx1, sx2, pn, ...
                gammapn, wz1gpn, wz2gpn, n, np);            
            for m=1:nterms
                for mp=1:nterms
                    if (mod(m+mp,2)==0)                        
                        Ky = sheatKernCompute(sigmay, lengthY, sy1, sy2, qm, ...
                            gammaqm, wz1gqm, wz2gqm, m, mp);
                        pn2qm2 = pn(n)^2 + qm(m)^2;
                        pnp2qmp2 = pn(np)^2 + qm(mp)^2;
                        K = K + (Kx.*Ky)/(pn2qm2*pnp2qmp2);
                    end
                end
            end
        end
    end
end

sK = cK*K;
sK = real(sK);
K = (poissKern.sensitivity^2)*sK;
K = real(K);


end





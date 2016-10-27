function [K, erfcz1, erfcz2] = srbfhKernComputeErfc(sigmax, lengthX, s1, s2, w, gamma, n)

% SRBFHKERNCOMPUTE Compute an SRBFH kernel.
% FORMAT

%
% COPYRIGHT : Mauricio A. Alvarez, 2016

% KERN

%%%%% Code with Faddeeva_erfc
sinS1 = sin(w(n)*s1);
bTerm = sigmax*gamma(n)/2;
argz2 = (s2-lengthX)/sigmax;
z2 = argz2 + bTerm;
argz1 =  s2/sigmax;
z1 = argz1 + bTerm;
erfcz2 = Faddeeva_erfc(z2);
erfcz1 = Faddeeva_erfc(z1);
vecXp = (sigmax*sqrt(pi)/2)*imag(exp(bTerm.^2)*exp(gamma(n)*s2).*...
    (erfcz2 - erfcz1));
K = sinS1*vecXp';
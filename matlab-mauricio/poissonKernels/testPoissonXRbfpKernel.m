clc
close all
clear

x = -2;
xp = 1;
y = -0.8;
yp = 1.0;
nterms = 5;
sigmax = 0.1;
sigmax2 = sigmax^2;
sigmay = 0.1;
sigmay2 = sigmay^2;
a = 10;
b = 10;

C = 4/(a*b);

knum = 0;
for n=1:nterms
    for m = 1:nterms
        pn = (n*pi)/a;
        qm = (m*pi)/b;
        pnqm = pn^2 + qm^2;        
        kterm = (1/pnqm)*sin(pn*x)*sin(qm*y)*...
            integral(@(xi) sin(pn*xi).*exp(-((xi - xp).^2)./sigmax2),0,a, 'RelTol', 1e-10)*...
            integral(@(eta) sin(qm*eta).*exp(-((eta - yp).^2)./sigmay2),0,b, 'RelTol', 1e-10);
        knum = knum + kterm;
    end
end
knum = C*knum;
%%% Now using the kernel
kernPoisson = kernCreate([x y], 'poisson');
kernRbfp = kernCreate([x y], 'rbfp');
kernPoisson.nTerms = nterms;
kernPoisson.inverseWidthX = 2/sigmax2;
kernPoisson.inverseWidthY = 2/sigmay2;
kernPoisson.lengthX = a;
kernPoisson.lengthY = b;
kernRbfp.inverseWidthX = 2/sigmax2;
kernRbfp.inverseWidthY = 2/sigmay2;
kanal = poissonXrbfpKernCompute(kernPoisson, kernRbfp, [x y], [xp yp]);





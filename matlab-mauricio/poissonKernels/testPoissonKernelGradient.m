clc
close all
clear

x = 1;
xp = 2;
y = 1.9;
yp = 1.0;
nterms = 5;
sigmax = 0.1;
sigmax2 = sigmax^2;
inverseWidthX = 2/sigmax2;
sigmay = 0.2;
sigmay2 = sigmay^2;
inverseWidthY = 2/sigmay2;
a = 10;
b = 10;

covGrad = 1;
epsilon = 1e-6;

kernPoisson = kernCreate([x y], 'poisson');
kernPoisson.nTerms = nterms;
kernPoisson.lengthX = a;
kernPoisson.lengthY = b;

% Numeric gradient with respect to inverseWidthX
kernPoisson.inverseWidthY = inverseWidthY;
kernPoisson.inverseWidthX = inverseWidthX + epsilon;
fplus = poissonKernCompute(kernPoisson, [x y], [xp yp]);
kernPoisson.inverseWidthX = inverseWidthX - epsilon;
fminus = poissonKernCompute(kernPoisson,[x y], [xp yp]);
gnum = 0.5*(fplus - fminus)/epsilon;

% Analityc gradient with respect to inverseWidthX
kernPoisson.inverseWidthX = inverseWidthX;
ganal = poissonKernGradient(kernPoisson, [x y], [xp yp], covGrad);

% Numeric gradient with respect to inverseWidthY
kernPoisson.inverseWidthX = inverseWidthX;
kernPoisson.inverseWidthY = inverseWidthY + epsilon;
fplus = poissonKernCompute(kernPoisson, [x y], [xp yp]);
kernPoisson.inverseWidthY = inverseWidthY - epsilon;
fminus = poissonKernCompute(kernPoisson, [x y], [xp yp]);
gnum = 0.5*(fplus - fminus)/epsilon;

% Analityc gradient with respect to inverseWidthY
kernPoisson.inverseWidthY = inverseWidthY;
ganal = poissonKernGradient(kernPoisson, [x y], [xp yp], covGrad);

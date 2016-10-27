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
inverseWidthX = 2/sigmax2;
sigmay = 0.2;
sigmay2 = sigmay^2;
inverseWidthY = 2/sigmay2;
a = 10;
b = 10;

covGrad = 1;
epsilon = 1e-6;

kernPoisson = kernCreate([x y], 'poisson');
kernRbfp = kernCreate([x y], 'rbfp');
kernPoisson.nTerms = nterms;
kernPoisson.lengthX = a;
kernPoisson.lengthY = b;

% Numeric gradient with respect to inverseWidthX
kernPoisson.inverseWidthY = inverseWidthY;
kernRbfp.inverseWidthY = inverseWidthY;
kernPoisson.inverseWidthX = inverseWidthX + epsilon;
kernRbfp.inverseWidthX = inverseWidthX + epsilon;
fplus = poissonXrbfpKernCompute(kernPoisson, kernRbfp, [x y], [xp yp]);
kernPoisson.inverseWidthX = inverseWidthX - epsilon;
kernRbfp.inverseWidthX = inverseWidthX - epsilon;
fminus = poissonXrbfpKernCompute(kernPoisson, kernRbfp, [x y], [xp yp]);
gnum = 0.5*(fplus - fminus)/epsilon;

% Analityc gradient with respect to inverseWidthX
kernPoisson.inverseWidthX = inverseWidthX;
kernRbfp.inverseWidthX = inverseWidthX;
ganal = poissonXrbfpKernGradient(kernPoisson, kernRbfp, [x y], [xp yp], covGrad);

% Numeric gradient with respect to inverseWidthY
kernPoisson.inverseWidthX = inverseWidthX;
kernRbfp.inverseWidthX = inverseWidthX;
kernPoisson.inverseWidthY = inverseWidthY + epsilon;
kernRbfp.inverseWidthY = inverseWidthY + epsilon;
fplus = poissonXrbfpKernCompute(kernPoisson, kernRbfp, [x y], [xp yp]);
kernPoisson.inverseWidthY = inverseWidthY - epsilon;
kernRbfp.inverseWidthY = inverseWidthY - epsilon;
fminus = poissonXrbfpKernCompute(kernPoisson, kernRbfp, [x y], [xp yp]);
gnum = 0.5*(fplus - fminus)/epsilon;

% Analityc gradient with respect to inverseWidthY
kernPoisson.inverseWidthY = inverseWidthY;
kernRbfp.inverseWidthY = inverseWidthY;
ganal = poissonXrbfpKernGradient(kernPoisson, kernRbfp, [x y], [xp yp], covGrad);

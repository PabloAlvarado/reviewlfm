clc
clear
close all

s = RandStream('mt19937ar', 'Seed', 1e2);
RandStream.setGlobalStream(s);
load('dataPoissonExample2.mat')

X1 = [xx1(:) xx2(:)];
lengthX = max(X1(:,1)) - min(X1(:,1));
lengthY = max(X1(:,2)) - min(X1(:,2));
y1 = ff_p(:);

options = multigpOptions('ftc');
options.optimiser = 'scg';
options.kernType = 'poisson';
options.kern.nTerms = 5;
options.kern.lenghtX = lengthX;
options.kern.lenghtY = lengthY;
options.nlf = 1;

X = cell(options.nlf+1,1);
y = cell(options.nlf+1,1);
XTest = cell(options.nlf+1,1);

X{1} = zeros(1,2);
y{1} = [];
options.bias(1) = 0;
options.scale(1) = 0;
X{2} = X1;
y{2} = y1;   
XTest{1} = X1;
XTest{2} = X1;    
options.bias(2) = mean(y1);
options.scale(2) = std(y1);


% Set the input and ouput dimensions
q = 2;
d = 1 + options.nlf;

warning('off','multiKernParamInit:noCrossKernel');
warning('off','multigp:FixedParam');

% Creates model
model = multigpCreate(q, d, X, y, options);

% Optimizes model
display = 1;
iters = 100;
trainingTime = cputime;
model = multigpOptimise(model, display, iters);
trainingTime = cputime - trainingTime;

[mu, varsigma] = multigpPosteriorMeanVar(model,  XTest);

figure
pcolor(xx1, xx2, reshape(ff_p, size(xx1)))
figure
pcolor(xx1, xx2, reshape(mu{2}, size(xx1)))
figure
surf(xx1, xx2, reshape(ff_p, size(xx1)))
figure
surf(xx1, xx2, reshape(mu{2}, size(xx1)))



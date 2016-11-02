clc
clear
close all
s = RandStream('mt19937ar','Seed',1e4);
RandStream.setGlobalStream(s);

N = 100;
nTrain = 10;
x = linspace(0,1, N)';
nSamples = 1;
nSamples2 = 3;
lengthScale = 0.2;
varWhite = 0.01;
kernType = 'rbf';
kernType2 = {'cmpnd', 'rbf', 'white'};
iters = 100;
display = true;

kern = kernCreate(x, 'rbf');
params = kernExtractParam(kern);
params(1) = log(1/(lengthScale^2));
kern = kernExpandParam(kern, params);
K = kernCompute(kern, x);
f = (gsamp(zeros(N,1), K, nSamples))';
% Adding noise 
y = f + 0.1*randn(N,1); 

% ds1 = rand('twister');
% ds2 = randn('seed');
index = randperm(N);
indexTrain = sort(index(1:nTrain));
%indexTrain = indexTrain([1 2 4 3]);
xTrain = x(indexTrain);
fTrain = f(indexTrain);
yTrain = y(indexTrain);
options = gpOptions('ftc');
options.kern = kernType2;
options.scale2var1 = true;
gpModel = gpCreate(1,1,xTrain, yTrain, options);
gpModel = gpOptimise(gpModel, display, iters);
% gpModel.kern.comp{1}.inverseWidth = 1/(lengthScale^2);
% gpModel.kern.comp{1}.variance = 1;
xTest = linspace(0, 1, N)';
%[mu, varsigma] = gpPosteriorMeanVar(gpModel, xTest);
[mu, covar] = gpPosteriorMeanCovar(gpModel, xTest);
varsigma = diag(covar{1});
fp = [(mu+2*real(sqrt(varsigma)));flip((mu-2*real(sqrt(varsigma))),1)];
a = fill([xTest; flip(xTest,1)], fp, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
hold on
b = plot(xTrain, yTrain, 'k.');
set(b, 'markersize', 25)
hold on
c = plot(xTest, mu,'k-');
set(c, 'lineWidth', 1.5);
d = plot(xTest, f, 'k--');
set(d, 'lineWidth', 1.5);
%%%
RandStream.setGlobalStream(s);
fPost = (gsamp(mu, covar{1}, nSamples2))';
% for i=1:nSamples2
%     hold on
%     d = plot(xTest, fPost(:,i),'k--');
%     set(d, 'lineWidth', 1);    
% end
%%%
xlabel('$t$', 'interpreter', 'latex', 'fontsize', 18)
ylabel('$f(t)$', 'interpreter', 'latex', 'fontsize', 18)
set(gca, 'fontname', 'arial', 'fontsize', 15, 'ylim', ylim, 'Color', 'none')
print('-depsc', '../pics/toyGPoutput', '-loose');



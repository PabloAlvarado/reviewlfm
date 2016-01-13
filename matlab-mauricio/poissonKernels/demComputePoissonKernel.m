clc
clear
close all

s = RandStream('mt19937ar', 'Seed', 1e2);
RandStream.setGlobalStream(s);

nX = 50;
nY = 50;
nterms = 5;
x = linspace(-1,1, nX)';
y = linspace(-1,1, nY)';
[X1, Y1] = meshgrid(x, y);
X = [X1(:) Y1(:)];
lengthX = max(x) - min(x);
lengthY = max(y) - min(y);

kern = kernCreate(X, 'poisson');
kern.nTerms = nterms;
kern.lengthX = lengthX;
kern.lengthY = lengthY;

% sigmax = 1;
% sigmay = 5;
% 
% kern.inverseWidthX = 2/sigmax^2;
% kern.inverseWidthY = 2/sigmay^2;


kern.inverseWidthX = 100;
kern.inverseWidthY = 10;
kern.sensitivity = 10*randn;

itime = cputime;
K = kernCompute(kern, X);
tpoint_true = cputime - itime;
imagesc(K)
diagK = kernDiagCompute(kern, X);
K = K + 1e-10*eye(size(K,1)); % Add a small regularization term
y = real(gsamp(zeros(size(K, 1),1), K, 1))'; % Get a sample from the joint GP
surf(X1, Y1, reshape(y, size(X1)))

% ntermsv = [1 3 5 10 15 20 30 50];
% Km = cell(length(ntermsv), 1);
% ym = cell(length(ntermsv), 1);
% tm = zeros(length(ntermsv), 1);
% 
% for n=1:length(ntermsv)
%     kern.nTerms = ntermsv(n);
%     itime = cputime;
%     Km{n} = kernCompute(kern, X);    
%     tm(n) = cputime - itime;
%     Kml = Km{n};
%     Kml = Kml + 1e-10*eye(size(Kml,1));
%     s = RandStream('mt19937ar', 'Seed', 1e8);
%     RandStream.setGlobalStream(s);
%     ym{n} = real(gsamp(zeros(size(Kml, 1),1), Kml, 1))';    
% end

% epsilon = 1e-6;
% 
% % Test the derivative for sigmax
% % param = sigmax + epsilon;
% 
% X = [0.2 -0.7];
% 
% % Test the derivative for inverseWidthX
% param = kern.inverseWidthX;
% kern.inverseWidthX = param + epsilon;
% K1 = kernCompute(kern, X);
% f1 = sum(sum(K1));
% %f1 = logdet(K1);
% kern.inverseWidthX = param - epsilon;
% K2 = kernCompute(kern, X);
% f2 = sum(sum(K2));
% %f2 = logdet(K2);
% gnum = 0.5*(f1-f2)/epsilon;
% kern.inverseWidthX = param;
% K = kernCompute(kern, X);
% %covGrad = (pdinv(K))';
% covGrad = ones(size(K));
% gana = kernGradient(kern, X, covGrad);

% Test the derivative for the sensitivity
% covGrad = ones(size(K));
% param = kern.sensitivity;
% kern.sensitivity = param + epsilon;
% K1 = kernCompute(kern, X);
% kern.sensitivity = param - epsilon;
% K2 = kernCompute(kern, X);
% gnum = 0.5*(K1-K2)/epsilon;
% kern.sensitivity = param;
% gana = kernGradient(kern, X, covGrad);







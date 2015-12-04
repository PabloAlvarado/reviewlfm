% DEMSAMPLESECONDORDERLFM

% LFM

clc
clear
close all
s = RandStream('mt19937ar', 'Seed', 1e8);
RandStream.setDefaultStream(s);

% Generate Samples from the model
N = 100; % Number of time points
t = linspace(0, 1, N)'; 
D = 1; % Number of outputs
%warning('off', 'multiKernParamInit:noCrossKernel') % Removes warning
kern = kernCreate(t, {'multi', 'rbf', 'lfm'});
lengthscale = 0.1;
massc = 1; % Mass constant
damperc = 40; % Damper constant = lambda if mass = 1
springc = 10; % Spring constant = gamma  if mass = 1
sensitivityc = 100; % Sensitivity parameter
% Initialize the lengthscale for the rbf kernel
kern.comp{1}.inverseWidth = 1/(lengthscale.^2);
kern.comp{1}.variance = 1;
for d=1:D
   kern.comp{1+d}.inverseWidth = 1/(lengthscale.^2);
   kern.comp{1+d}.damper = damperc(d);
   kern.comp{1+d}.spring = springc(d);
   kern.comp{1+d}.mass = massc(d);
   kern.comp{1+d}.sensitivity = sensitivityc(d);
end
K = kernCompute(kern, t); % Compute the kernel
K = K + 1e-10*eye(size(K,1)); % Add a small regularization term
yV = real(gsamp(zeros(size(K, 1),1), K, 1))'; % Get a sample from the joint GP
% Get the data and plot it
f = yV(N+1:end);
u = yV(1:N);
% First we plot the output
figure
plot(t, f, 'linewidth', 1.5)
xlabel('Time')
ylabel('Output function')
title('Sample from the output function')
% Now we plot the input
figure
plot(t, u, 'linewidth', 1.5)
xlabel('Time')
ylabel('Latent function')
title('Sample for the input function')
% Now, we plot the covariance matrix Kff
figure
imagesc(t, t, K(N+1:end, N+1:end))
title('Covariance matrix 2nd order LFM one output')
% We plot the covariance for the latent function, Kuu 
figure
imagesc(t, t, K(1:N, 1:N))
title('Covariance matrix latent function or latent force')
% Now, we plot the covariance between the output and the latent force, Kfu
figure
imagesc(t, t, K(N+1:end, 1:N))
title('Cross Covariance matrix between output and the latent force')

% Given some measurements with noise, we now infer the force
Nk = 30;
index = randperm(N);
indexk = sort(index(1:Nk));
fk = f(indexk);
tk = t(indexk);
% We add some noise to the measurements fk,
sigma2 = 0.01*var(f);
yk = fk +  sqrt(sigma2)*randn(Nk, 1);
close all
plot(tk, yk, '.r', t, f, 'b', 'markersize', 15, 'linewidth', 1.5 )

% We now use the measurements to compute the posterior over the latent
% force
Kaux = kernCompute(kern, tk);
Kfkfk = Kaux(Nk+1:end, Nk+1:end);
Kykyk = Kfkfk + sigma2*eye(Nk);
invKykyk = Kykyk\(eye(Nk));
Kaux = kernCompute(kern, t, tk);
Kustar_f = Kaux(1:N, Nk+1:end);
Kustar_ustar = kernCompute(kern.comp{1}, t);
mean_pred_u = Kustar_f*invKykyk*yk;
var_pred_u = diag(Kustar_ustar - Kustar_f*invKykyk*(Kustar_f'));

figure
g = [(mean_pred_u + 2*real(sqrt(var_pred_u))); flipdim((mean_pred_u-2*real(sqrt(var_pred_u))),1)];
fill([t; flipdim(t,1)], g, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
hold on
plot(t, mean_pred_u, 'k--', 'linewidth', 1.5)
hold on
plot(t, u, 'k','linewidth', 1.5)
set(gca, 'Color', 'none')

save('dataODEwithRBF.mat')

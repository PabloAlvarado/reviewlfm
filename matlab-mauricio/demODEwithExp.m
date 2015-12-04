% DEMSAMPLESECONDORDERLFM

% LFM


%magnSigma2 = 1;
%lengthScale = 1;
%matern = @(tau) magnSigma2 * ...
%    exp(-sqrt(3)*abs(tau)./lengthScale).*(sqrt(3)*abs(tau)/lengthScale+1);

clc
clear
close all
s = RandStream('mt19937ar', 'Seed', 1e4);
RandStream.setDefaultStream(s);

% Generate Samples from the model
N = 100; % Number of time points
t = linspace(0, 1, N)'; 
D = 1; % Number of outputs
%warning('off', 'multiKernParamInit:noCrossKernel') % Removes warning
kern = kernCreate(t, {'multi', 'rbf', 'lfm'});
lengthscale = 0.1;
inverseWidth = 1/lengthscale;
mass = 1; % Mass constant
damper = 40; % Damper constant = lambda if mass = 1
spring = 10; % Spring constant = gamma  if mass = 1
sensitivity = 100; % Sensitivity parameter
% Compute the kernel Kuu
Kuu = computeExpKernel(t, t, inverseWidth);
Kff = computeLfmExpKernel(t, t, mass, damper, spring, sensitivity, inverseWidth);
Kfu = computeLfmXExpKernel(t, t, mass, damper, spring, sensitivity, inverseWidth);

K = [Kuu Kfu'; Kfu Kff];

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
sigma2 = 0.00000001*var(f);
yk = fk +  sqrt(sigma2)*randn(Nk, 1);
close all
plot(tk, yk, '.r', t, f, 'b', 'markersize', 15, 'linewidth', 1.5 )

% We now use the measurements to compute the posterior over the latent
% force
Kfkfk = computeLfmExpKernel(tk, tk, mass, damper, spring, sensitivity, inverseWidth);
Kykyk = Kfkfk + sigma2*eye(Nk);
invKykyk = Kykyk\(eye(Nk));
Kustar_f = (computeLfmXExpKernel(tk, t, mass, damper, spring, sensitivity, inverseWidth))';
Kustar_ustar = computeExpKernel(t, t, inverseWidth);
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
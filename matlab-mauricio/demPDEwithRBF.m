% DEMSAMPLESECONDORDERLFM

% LFM

clc
clear
close all
s = RandStream('mt19937ar', 'Seed', 1e8);
RandStream.setDefaultStream(s);

nterms = 5;  % Number of terms in the Fourier solutions

if exist(['dataPDEwithRBF' num2str(nterms) 'terms.mat'], 'file')
    load(['dataPDEwithRBF' num2str(nterms) 'terms.mat'])
else
    % Generate Samples from the model
    Nx = 40; % Number of space points
    Nt = 40; % Number of time points
    t = linspace(0, 1, Nt)';
    lengthX = 5; % Length of domain in X
    x = linspace(0, lengthX, Nx)';
    % We now organize the input data: first column time dimension, second
    % column space dimension
    startVal = 1;
    endVal = 0;
    tx = zeros(Nt*Nx, 2);
    for i=1:Nt
        endVal = endVal + Nx;
        aux = [t(i)*ones(Nx,1) x];
        tx(startVal:endVal, 1:2) = aux;
        startVal = endVal + 1;
    end
    %warning('off', 'multiKernParamInit:noCrossKernel') % Removes warning
    kern = kernCreate(tx, {'multi', 'rbfh', 'heat'});
    lengthscale_space = 0.5;
    lengthscale_time = 0.1;
    decay = 2;  % Parameters \lambda in the equation
    diffusion = 1e-3; % Diffusion constant
    sensitivity = 1; % Sensitivity parameter
    % Initialize the lengthscales for the rbf kernels
    kern.comp{1}.inverseWidthSpace = 1/(lengthscale_space.^2);
    kern.comp{1}.inverseWidthTime = 1/(lengthscale_time.^2);
    kern.comp{2}.inverseWidthSpace = 1/(lengthscale_space.^2);
    kern.comp{2}.inverseWidthTime = 1/(lengthscale_time.^2);
    kern.comp{2}.decay = decay;
    kern.comp{2}.diffusion = diffusion;
    kern.comp{2}.sensitivity = sensitivity;
    kern.comp{2}.lengthX = lengthX;
    K = kernCompute(kern, tx); % Compute the kernel
    K = K + 1e-10*eye(size(K,1)); % Add a small regularization term
    yV = real(gsamp(zeros(size(K, 1),1), K, 1))'; % Get a sample from the joint GP
    % We separate the kernels
    Kff = K(Nx*Nt+1:end, Nx*Nt+1:end);
    Kuu = K(1:Nx*Nt, 1:Nx*Nt);
    Kuf = K(1:Nx*Nt, Nx*Nt+1:end);
end
% Now, we plot the covariance matrix Kff
figure
imagesc(t, x, Kff)
title('Covariance matrix PDE one output')
% We plot the covariance for the latent function, Kuu
figure
imagesc(t, x, Kuu)
title('Covariance matrix latent function or latent force')
% Now, we plot the covariance between the output and the latent force, Kfu
figure
imagesc(t, x, Kuf)
title('Cross Covariance matrix between output and the latent force')
% Get the data and plot it
f = yV(Nt*Nx+1:end);
u = yV(1:Nt*Nx);
[tgrid, xgrid] = meshgrid(t, x);
F = reshape(f, size(tgrid));
U = reshape(u, size(tgrid));
% First we plot the output
figure
surf(tgrid, xgrid, F)
xlabel('Time')
ylabel('Space')
zlabel('Output function')
title('Sample from the output function')
% Now we plot the input
figure
surf(tgrid, xgrid, U)
xlabel('Time')
ylabel('Space')
zlabel('input function')
title('Sample for the input function')

% Given some measurements with noise, we now infer the force
percent = 0.5; % Percentage of data used to infer the latent force
Nk = floor(0.5*Nt*Nx);
index = randperm(Nt*Nx);
indexk = sort(index(1:Nk));
fk = f(indexk);
txk = tx(indexk, :);
% We add some noise to the measurements fk,
sigma2 = 0.01*var(f);
yk = fk +  sqrt(sigma2)*randn(Nk, 1);
close all
surf(tgrid, xgrid, F)
hold on
plot3(txk(:,1), txk(:, 2), yk,  '.k', 'markersize', 15, 'linewidth', 1.5 )

% We now use the measurements to compute the posterior over the latent
% force
Kaux = kernCompute(kern, txk);
Kfkfk = Kaux(Nk+1:end, Nk+1:end);
Kykyk = Kfkfk + sigma2*eye(Nk);
invKykyk = Kykyk\(eye(Nk));
Kaux = heatXrbfhKernCompute(kern.comp{2}, kern.comp{1}, txk, tx);
Kustar_f = Kaux';
Kustar_ustar = kernCompute(kern.comp{1}, tx);
mean_pred_u = Kustar_f*invKykyk*yk;
var_pred_u = diag(Kustar_ustar - Kustar_f*invKykyk*(Kustar_f'));
% Plot of the real latent force
figure
contourf(tgrid, xgrid, U)
colorbar
title('Real latent force')
% Plot of the predicted mean of the latent force
figure
contourf(tgrid, xgrid, reshape(mean_pred_u, size(tgrid)))
colorbar
title('Predicted mean of the latent force')
% Plot of the predicted variance of the latent force
figure
contourf(tgrid, xgrid, reshape(var_pred_u, size(tgrid)))
colorbar
title('Predicted variance of the latent force')

mse_error = mean((mean_pred_u - u).^2);

save(['dataPDEwithRBF' num2str(nterms) 'terms.mat'])

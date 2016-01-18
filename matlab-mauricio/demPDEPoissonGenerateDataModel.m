% DEMSAMPLESECONDORDERLFM

% LFM

clc
clear
close all
s = RandStream('mt19937ar', 'Seed', 1e8);
RandStream.setGlobalStream(s);

nterms = 10;  % Number of terms in the Fourier solutions

if exist(['dataPoissonExample2' num2str(nterms) 'terms.mat'], 'file')
    load(['dataPoissonExample2' num2str(nterms) 'terms.mat'])
else
    % Generate Samples from the model
    Nx = 30; % Number of space points in X
    Ny = 30; % Number of space points in Y
    lengthX = 2; % Length of domain in X
    x = linspace(-1, 1, Nx)';
    lengthY = 2; % Length of domain in X
    y = linspace(-1, 1, Ny)';
    % We now organize the input data
    [X1, Y1] = meshgrid(x, y);
    xy = [X1(:) Y1(:)];    
    %warning('off', 'multiKernParamInit:noCrossKernel') % Removes warning
    kern = kernCreate(xy, {'multi', 'rbfp', 'poisson'});
    lengthscale_spaceX = 0.1;
    lengthscale_spaceY = 0.5;    
    sensitivity = 1; % Sensitivity parameter
    % Initialize the lengthscales for the rbf kernels
    kern.comp{1}.inverseWidthX = 1/(lengthscale_spaceX.^2);
    kern.comp{1}.inverseWidthY = 1/(lengthscale_spaceY.^2);
    kern.comp{2}.inverseWidthX = 1/(lengthscale_spaceX.^2);
    kern.comp{2}.inverseWidthY = 1/(lengthscale_spaceY.^2);        
    kern.comp{2}.sensitivity = sensitivity;
    kern.comp{2}.lengthX = lengthX;
    kern.comp{2}.lengthY = lengthY;
    kern.comp{2}.nTerms = nterms;
    K = kernCompute(kern, xy); % Compute the kernel
    K = K + 1e-10*eye(size(K,1)); % Add a small regularization term
    yV = real(gsamp(zeros(size(K, 1),1), K, 1))'; % Get a sample from the joint GP
    % We separate the kernels
    Kff = K(Nx*Ny+1:end, Nx*Ny+1:end);
    Kuu = K(1:Nx*Ny, 1:Nx*Ny);
    Kuf = K(1:Nx*Ny, Nx*Ny+1:end);
end
% Now, we plot the covariance matrix Kff
figure
pcolor(Kff)
title('Covariance matrix PDE Poisson one output')
% We plot the covariance for the latent function, Kuu
figure
pcolor(Kuu)
title('Covariance matrix latent function or latent force')
% Now, we plot the covariance between the output and the latent force, Kfu
figure
pcolor(Kuf)
title('Cross Covariance matrix between output and the latent force')
% Get the data and plot it
f = yV(Nx*Ny+1:end);
u = yV(1:Nx*Ny);
% First we plot the output
figure
surf(X1, Y1, reshape(f, size(X1)))
xlabel('Space X')
ylabel('Space Y')
zlabel('Output function')
title('Sample from the output function')
% Now we plot the input
figure
surf(X1, Y1, reshape(u, size(X1)))
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

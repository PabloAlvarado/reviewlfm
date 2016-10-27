% DEMSAMPLESECONDORDERLFM

% LFM

clc
clear
close all
s = RandStream('mt19937ar', 'Seed', 1e5);
RandStream.setGlobalStream(s);

nterms = 5;  % Number of terms in the Fourier solutions

% if exist(['dataPoissonExample2' num2str(nterms) 'terms.mat'], 'file')
%     load(['dataPoissonExample2' num2str(nterms) 'terms.mat'])
% else
    % Generate Samples from the model
    Nx = 30; % Number of space points in X
    Ny = 30; % Number of space points in Y
    lengthX = 1; % Length of domain in X
    x = linspace(0, 1, Nx)';
    lengthY = 1; % Length of domain in X
    y = linspace(0, 1, Ny)';
    % We now organize the input data
    [X1, Y1] = meshgrid(x, y);
    xy = [X1(:) Y1(:)];    
    %warning('off', 'multiKernParamInit:noCrossKernel') % Removes warning
    kern = kernCreate(xy, {'multi', 'rbfp', 'poisson'});
    lengthscale_spaceX = 0.3;
    lengthscale_spaceY = 0.3;    
    sensitivity = 1; % Sensitivity parameter
    % Initialize the lengthscales for the rbf kernels
    kern.comp{1}.inverseWidthX = 2/(lengthscale_spaceX.^2);
    kern.comp{1}.inverseWidthY = 2/(lengthscale_spaceY.^2);
    kern.comp{2}.inverseWidthX = 2/(lengthscale_spaceX.^2);
    kern.comp{2}.inverseWidthY = 2/(lengthscale_spaceY.^2);        
    kern.comp{2}.sensitivity = sensitivity;
    kern.comp{2}.lengthX = lengthX;
    kern.comp{2}.lengthY = lengthY;
    kern.comp{2}.nTerms = nterms;
    K = kernCompute(kern, xy); % Compute the kernel
    %K = K + 1e-5*eye(size(K,1)); % Add a small regularization term
    yV = real(gsamp(zeros(size(K, 1),1), K, 1))'; % Get a sample from the joint GP
    % We separate the kernels
    Kff = K(Nx*Ny+1:end, Nx*Ny+1:end);
    Kuu = K(1:Nx*Ny, 1:Nx*Ny);
    Kuf = K(1:Nx*Ny, Nx*Ny+1:end);
%end
% Now, we plot the covariance matrix Kff
figure
imagesc(Kff)
title('Covariance matrix PDE Poisson one output')
% We plot the covariance for the latent function, Kuu
figure
imagesc(Kuu)
title('Covariance matrix latent function or latent force')
% Now, we plot the covariance between the output and the latent force, Kfu
figure
imagesc(Kuf)
title('Cross Covariance matrix between output and the latent force')
% Get the data and plot it
f = yV(Nx*Ny+1:end);
u = yV(1:Nx*Ny);
F = reshape(f, size(X1));
U = reshape(u, size(X1));
% First we plot the output
figure
surf(X1, Y1, F)
xlabel('Space X')
ylabel('Space Y')
zlabel('Output function')
title('Sample from the output function')
% Now we plot the input
figure
surf(X1, Y1, U)
xlabel('Time')
ylabel('Space')
zlabel('input function')
title('Sample for the input function')


% Given some measurements with noise, we now infer the force
percent = 0.5; % Percentage of data used to infer the latent force
Nk = floor(0.5*Nx*Ny);
index = randperm(Nx*Ny);
indexk = sort(index(1:Nk));
fk = f(indexk);
xyk = xy(indexk, :);
% We add some noise to the measurements fk,
%sigma2 = 0.001*var(f);
sigma2= 0;
yk = fk +  sqrt(sigma2)*randn(Nk, 1);
%close all
figure
surf(X1, Y1, reshape(f, size(X1)))
hold on
plot3(xyk(:,1), xyk(:, 2), yk,  '.k', 'markersize', 15, 'linewidth', 1.5 )

% savedata = true;
% xx1 = X1;
% xx2 = Y1;
% uu = u;
% ff_p = f;
% %X1 = xyk;
% if savedata
%     save('./dataPoissonExample2.mat', 'xx1', 'xx2', 'uu', 'ff_p')
% end

% % We now use the measurements to compute the posterior over the latent
% % force
Kaux = kernCompute(kern, xyk);
Kfkfk = Kaux(Nk+1:end, Nk+1:end);
%Kykyk = Kfkfk + sigma2*eye(Nk);
Kykyk = Kfkfk + 1e-10*eye(Nk);
invKykyk = Kykyk\(eye(Nk));
Kaux = poissonXrbfpKernCompute(kern.comp{2}, kern.comp{1}, xyk, xy);
Kustar_f = Kaux';
Kustar_ustar = kernCompute(kern.comp{1}, xy);
mean_pred_u = Kustar_f*invKykyk*yk;
var_pred_u = diag(Kustar_ustar - Kustar_f*invKykyk*(Kustar_f'));
% Plot of the real latent force
figure
contourf(X1, Y1, U)
colorbar
title('Real latent force')
figure
surf(X1, Y1, reshape(mean_pred_u, size(X1)))
colorbar
title('Predicted mean of the latent force')
% Plot of the predicted mean of the latent force
figure
contourf(X1, Y1, reshape(mean_pred_u, size(X1)))
colorbar
title('Predicted mean of the latent force')
% Plot of the predicted variance of the latent force
figure
contourf(X1, Y1, reshape(var_pred_u, size(X1)))
colorbar
title('Predicted variance of the latent force')

mse_error = mean((mean_pred_u - u).^2);

% save(['dataPDEwithRBF' num2str(nterms) 'terms.mat'])

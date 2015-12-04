% DEMSAMPLESECONDORDERLFM

% LFM

clc
clear
close all
s = RandStream('mt19937ar', 'Seed', 1e8);
RandStream.setDefaultStream(s);

nterms = 5;  % Number of terms in the Fourier solutions, choose 5 or 10

load(['dataPDEwithRBF' num2str(nterms) 'terms_NoK.mat'])

close all

% Plot the data and the true output
surf(tgrid, xgrid, F)
hold on
plot3(txk(:,1), txk(:, 2), yk,  '.k', 'markersize', 15, 'linewidth', 1.5 )

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


%
% Simulate data from the model using a finite differences
% approximation
%

    %%
    % Define the covariance function and draw a sample
    %
    lengthscale_time = 0.1;
    lengthscale_space = 0.5;
  
    se_2d = @(x,t) exp(-x.^2/lengthscale_space^2/2-t.^2/lengthscale_time^2/2);

    % Generate Samples from the model
    Nx = 40; % Number of space points
    Nt = 40; % Number of time points
    t = linspace(0, 1, Nt)';
    lengthX = 5; % Length of domain in X
    x = linspace(0, lengthX, Nx)';

    [tgrid, xgrid] = meshgrid(t, x);

    K = zeros(length(x));
    for i=1:length(xgrid(:))
        for j=1:length(xgrid(:))
            K(i,j) = se_2d(xgrid(i) - xgrid(j), tgrid(i) - tgrid(j));
        end
    end
    
    [evec, eval] = eig(K);

    deig = diag(eval);

    if (~isreal(deig)) | any(deig<0), 
        warning('Covariance Matrix is not OK, redefined to be positive definite');
        eval=abs(eval);
    end

    coeffs = randn(1, size(K,1))*sqrt(eval);

    u = coeffs*evec';
    u = u(:);
    
    clf;
    pcolor(tgrid,xgrid,reshape(u, size(xgrid)));
    shading flat
    colorbar;
    
    %%
    % Solve the PDE with finite-differences together
    % with a backward Euler
    %
    decay = 2;  % Parameters \lambda in the equation
    diffusion = 1e-3; % Diffusion constant
    
    dt = t(2) - t(1);
    dx = x(2) - x(1);
    
    n = length(x);
    a = repmat(-2*diffusion/dx^2-decay,n,1);
    b = repmat(diffusion/dx^2,n,1);
    F = spdiags([b a b],-1:1,n,n);

    spy(F)

    uu = reshape(u,size(xgrid));
    f = zeros(size(xgrid));
    for k=2:length(t)
        % x_{k+1} = x_k + F x_{k+1} dt + u_{k+1} dt
        % x_{k+1} = (I - F dt)^-1 (x_k + u_{k+1} dt)
        f(:,k) = (eye(size(F,1)) - F * dt) \ (f(:,k-1) + uu(:,k) * dt);
    end
    
    f = f(:);
    
    clf;
    pcolor(tgrid,xgrid,reshape(f, size(xgrid)));
    colorbar;
    shading flat
    
    %%
    % Draw the measurements
    %    
    startVal = 1;
    endVal = 0;
    tx = zeros(Nt*Nx, 2);
    for i=1:Nt
        endVal = endVal + Nx;
        aux = [t(i)*ones(Nx,1) x];
        tx(startVal:endVal, 1:2) = aux;
        startVal = endVal + 1;
    end
    
    Nk = floor(0.5*Nt*Nx);
    
    index = randperm(Nt*Nx);
    indexk = sort(index(1:Nk));
    fk = f(indexk);
    txk = tx(indexk, :);
    sigma2 = 0.01*var(f);
    yk = fk + sqrt(sigma2)*randn(Nk, 1);

    FF = reshape(f, size(tgrid));
    UU = reshape(u, size(tgrid));
    
    % First we plot the output
    clf;
    subplot(1,2,1);
    surf(tgrid, xgrid, FF)
    hold on
    plot3(txk(:,1), txk(:, 2), yk,  '.k', 'markersize', 15, 'linewidth', 1.5 )
    xlabel('Time')
    ylabel('Space')
    zlabel('Output function')
    title('Sample from the output function')

    % Now we plot the input
    subplot(1,2,2);
    surf(tgrid, xgrid, UU)    
    xlabel('Time')
    ylabel('Space')
    zlabel('input function')
    title('Sample for the input function')

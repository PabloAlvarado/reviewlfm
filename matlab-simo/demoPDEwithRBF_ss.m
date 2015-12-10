%
% State-space (Kalman) PDE solution with SE covariance
%

    %%
    % Load the data
    %
%    load ../matlab-mauricio/dataPDEwithRBF10terms_NoK.mat
%    load ../matlab-mauricio/dataPDEwithRBF5terms_NoK.mat
%    load dataPDEwithRBF.mat

%    load ../matlab-mauricio/dataPDEwithRBF.mat
    load ../matlab-mauricio/dataRBFwithPDE_SIMFD.mat
    mean_pred_u = mean_pred_u10;
    
    %%
    % Plot output
    %
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

    %%
    % Temporal part in state-space form
    %
%    s_t = 1; % Scale
%    ell_t = 0.1;    
    s_t = 1
    ell_t = lengthscale_time
%    order_t = 4; % Order of state-space model
    
    se_cov_t = @(t) s_t^2 * exp(-t.^2/2/ell_t^2)

    [Bte,Ate] = se_pade(4,8,s_t,ell_t)
%    [Bte,Ate] = se_taylor(order_t,s_t,ell_t);
    [Fte,Lte,qte,Hte] = ratspec_to_ss(Bte,Ate);
    [Fte,Lte,Hte] = ss_balance(Fte,Lte,Hte);

    Fte
    Lte
    Hte
    qte
    
    tau = -5*ell_t:ell_t/10:5*ell_t;
    c = ss_cov(tau,Fte,Lte,qte,Hte);
    cond(Fte)

    Pte = lyapchol(Fte,Lte*sqrt(qte));
    Pte = Pte' * Pte;
    
    clf;
    plot(tau,se_cov_t(tau),tau,c,'r--');
    title('Accuracy of the state-space approximation of SE');
    legend('Exact SE','State-space Approx');
    grid on;    
    
    %%
    % Form a series approximation to the spatial covariance
    %
    
    % Spatial covariance function
%    s_x = 1; % Scale
%    ell_x = 0.1;    
    s_x = 1
    ell_x = lengthscale_space
    
    se_cov_x = @(x,y) s_x^2 * exp(-(x-y).^2/2/ell_x^2)
    se_spec_x = @(w) s_x^2 * sqrt(2*pi) * ell_x * exp(-ell_x^2 * w.^2/2);
    
    N = 20; % How many series terms
    LL = lengthX;  % Size of the domain

     % Laplace operator eigenvalues and functions for [-LL,LL]
%
%    eigenval = @(n) (n(:)'*pi/2/LL).^2;
%    eigenfun = @(n,x) LL^(-1/2)*sin(kron(n(:)'*pi,(x(:)+LL)/2/LL));

    % Laplace operator eigenvalues and functions for [0,L]
    eigenval = @(n) (n(:)'*pi/LL).^2
    eigenfun = @(n,x) (LL/2)^(-1/2)*sin(kron(n(:)'*pi,x(:)/LL))
    eigenfun_int = @(n,x) (LL/2)^(-1/2)*sin(n*pi*x/LL)

    % Project the spectral density onto the basis
    %{
    Qc = zeros(N);
    for i=1:N
        fprintf('Computing noise projections %d/%d\n',i,N);
        for j=1:N
            Qc(i,j) = integral2(@(x,y) eigenfun_int(i,x) .* se_cov_x(x,y) .* eigenfun_int(j,y),0,LL,0,LL);
        end
    end
    %}
    
    % This is a bit rougher approximation
    se_spec_x(sqrt(eigenval(3)))
    Qc = zeros(N);
    for i=1:N
        Qc(i,i) = se_spec_x(sqrt(eigenval(i)));
    end

    
    %%
    % Show the approximate covariance (not needed in estimation)
    %
    %{
    xx = 0:0.1:LL;
    K = zeros(length(xx));
    Ka = zeros(length(xx));
    for i=1:size(K,1)
        fprintf('Computing approximate covariance %d/%d\n',i,size(K,1));
        for j=1:size(K,2)
            tx = xx(i);
            ty = xx(j);
            K(i,j) = se_cov_x(tx,ty);
            
            for ii=1:N
                for jj=1:N
                    Ka(i,j) = Ka(i,j) + Qc(ii,jj) * eigenfun(ii,tx) * eigenfun(jj,ty);
                end
            end
        end
    end

    clf;

    subplot(1,2,1);
    colormap(hot);
    [tmp1,tmp2] = meshgrid(xx);
    pcolor(tmp1,tmp2,K);
    shading flat
    colorbar;
    
    subplot(1,2,2);
    colormap(hot);
    [tmp1,tmp2] = meshgrid(xx);
    pcolor(tmp1,tmp2,Ka);
    shading flat
    colorbar;
    %}
    
    %%
    % Form the joint state-space model
    %

    % Evaluate the basis at measurements (for testing)
    ut = unique(txk(:,1));
    ind = find(txk(:,1) == ut(1));
    meas_x = txk(ind,2);
    
    E = zeros(length(meas_x),N);
    for i=1:length(meas_x)
        for j=1:N
            tx = meas_x(i);
            E(i,j) = eigenfun(j,tx);
        end
    end

    % Evaluate the basis at prediction points
    Ep = zeros(length(x),N);
    for i=1:length(x)
        for j=1:N
            tx = x(i);
            Ep(i,j) = eigenfun(j,tx);
        end
    end    
    
    % State-space model for the force
    Fgp = kron(eye(N),Fte);
    Lgp = sqrt(qte) * kron(eye(N),Lte);
    Hgpc = kron(eye(N),Hte); % Gives plain coefficients
    Hgp  = kron(E,Hte);      % Projects into measurement points
    Hgpp = kron(Ep,Hte);     % Projects into prediction points

    Pgp = lyap(Fgp,Lgp*Qc*Lgp');
    
    % State-space model for the PDE
    D = diffusion
    lam = decay
    
    Fpd = zeros(N);
    Lpd = eye(N);
    
    for i=1:N
        Fpd(i,i) = -lam + D * eigenval(i);
    end

    Hpdc = eye(N); % Gives plain coefficients
    Hpd  = E;       % Projects into measurement points
    Hpdp = Ep;      % Projects into prediction points
     
    % Form the joint state-space model
%    sigma2 = 0.1;
    
    F = blkdiag(Fpd,Fgp);
    d = size(Fpd,1);
    F(1:d,d+1:end) = Lpd * Hgpc;
    L = [zeros(d,size(Lgp,2)); Lgp];
    H = [Hpd zeros(size(Hpd,1),size(Fgp,1))];
    Hu = [zeros(size(Hgp,1),d) Hgp];
    Hp = [Hpdp zeros(size(Hpdp,1),size(Fgp,1))];
    Hup = [zeros(size(Hgpp,1),d) Hgpp];
    
%    [F,L,H,T] = ss_balance(F,L,H);
%    Hu = Hu * T;
    cond(F)
    
    clf;
    subplot(2,2,1);
    spy(F);
    
    subplot(2,2,2);
    spy(L);
    
    subplot(2,2,3);
    spy(H);

    subplot(2,2,4);
    spy(Hu);
    
    %%
    % Compute the discretization (we have a uniform
    % time grid so we can do this beforehand)
    %
    fprintf('Discretizing...\n');
    dt = t(2)-t(1);
    A = expm(F*dt);
    P = lyap(F,L*Qc*L');
    Q = P - A*P*A';

    fprintf('Done.\n');

    clf;
    subplot(2,2,1);
    spy(A);
    
    subplot(2,2,2);
    spy(Q);
    
    subplot(2,2,3);
    spy(Hu);

    subplot(2,2,4);
    spy(Hup);
    
    %%
    % Kalman filter
    %
    clf;
    kf_x_mu = zeros(length(x),length(t));
%    kf_x_V  = [];

    kf_u_mu = zeros(length(x),length(t));
%    kf_u_V  = [];
    
    m = zeros(size(F,1),1);
    %m(1) = x0;
    %m(2) = dx0;
    P = blkdiag(0*eye(size(Fpd,1)),Pgp);
    P0 = P;

    MM = zeros(size(m,1),length(t));
    PP = zeros(size(P,1),size(P,2),length(t));
    
    plot_list = length(t);
    
    tic;
    fprintf('Evaluating Hs...\n');
    for k=1:length(t)
        % Evaluate the basis at measurement locations
        % for each time
        ind = find(abs(txk(:,1) - t(k)) < eps);
        meas_x = txk(ind,2);
        y = yk(ind);
        
%        E = zeros(length(meas_x),N);
%        for i=1:length(meas_x)
%            xx = meas_x(i);
%            %for j=1:N
%            %    E(i,j) = eigenfun(j,xx);
%            %end
%        end
        E = eigenfun(1:N,meas_x);
        Hpd  = E;       % Projects into measurement points
        H = [Hpd zeros(size(Hpd,1),size(Fgp,1))];
        HH{k} = H;
    end
    toc
    
    fprintf('Running Kalman filter...\n');
    tic;
    
    for k=1:length(t)
        
        if k > 1
            m = A*m;
            P = A*P*A' + Q;
        end

        %{
        % Evaluate the basis at measurements
        ind = find(abs(txk(:,1) - t(k)) < eps);
        meas_x = txk(ind,2);
        y = yk(ind);
    
        if length(y) > 0
            E = zeros(length(meas_x),N);
            for i=1:length(meas_x)
                for j=1:N
                    xx = meas_x(i);
                    E(i,j) = eigenfun(j,xx);
                end
            end
            Hpd  = E;       % Projects into measurement points
            H = [Hpd zeros(size(Hpd,1),size(Fgp,1))];
        
            R = sigma2 * eye(length(y));
            S = H*P*H' + R;
            K = P * H' / S;
            m = m + K * (y - H*m);
            P = P - K * S * K';
            
        end
        %}
        ind = find(abs(txk(:,1) - t(k)) < eps);
        y = yk(ind);
        
        H = HH{k};
        R = sigma2 * eye(length(y));
        S = H*P*H' + R;
        K = P * H' / S;
        m = m + K * (y - H*m);
        P = P - K * S * K';
        
        MM(:,k) = m;
        PP(:,:,k) = P;
        
        kf_x_mu(:,k) = Hp*m;
%        kf_x_V  = [kf_x_V  Hp*P*Hp'];

        kf_u_mu(:,k) = Hup*m;
%        kf_u_V  = [kf_u_V  Hup*P*Hu'];

        if any(plot_list == k)
            subplot(2,2,1);
            pcolor(xgrid,tgrid,kf_x_mu);
            shading interp;
            title('KF f');
            
            subplot(2,2,2);
            pcolor(xgrid,tgrid,FF);
            shading interp;
            title('GP f');
            
            subplot(2,2,3);
            pcolor(xgrid,tgrid,kf_u_mu);
            shading interp;
            title('KF u');
            
            subplot(2,2,4);
            pcolor(xgrid,tgrid,UU);
            shading interp;
            title('GP u');

            drawnow;
        end

        %pause;
    end
    toc

    %%
    % Smoother
    %    
    MS = MM;
    PS = PP;
    ms = MM(:,end);
    Ps = PP(:,:,end);
    
%    ks_x_mu = zeros(length(x),length(t));
    ks_x_mu = kf_x_mu;
%    ks_x_V  = [];

%    ks_u_mu = zeros(length(x),length(t));
    ks_u_mu = kf_u_mu;
%    ks_u_V  = [];

    ks_x_mu(:,end) = Hp*ms;
    ks_u_mu(:,end) = Hup*ms;
    
    plot_list = 1;
    
    fprintf('Running smoother...\n');

    tic;
    for k=size(MM,2)-1:-1:1
        m = MM(:,k);
        P = PP(:,:,k);
        
        mp = A*m;
        Pp = A*P*A' + Q;
        
        G = P*A'/Pp;
        ms = m + G * (ms - mp);
        Ps = P + G * (Ps - Pp) * G';
        
        MS(:,k) = ms;
        PS(:,:,k) = Ps;
        
        ks_x_mu(:,k) = Hp*ms;
        ks_u_mu(:,k) = Hup*ms;
        
        if any(plot_list == k)
            subplot(2,2,1);
            pcolor(xgrid,tgrid,ks_x_mu);
            shading interp;
            title('KS f');
            
            subplot(2,2,2);
            pcolor(xgrid,tgrid,FF);
            shading interp;
            title('GP f');
            
            subplot(2,2,3);
            pcolor(xgrid,tgrid,ks_u_mu);
            shading interp;
            title('KS u');
            
            subplot(2,2,4);
            pcolor(xgrid,tgrid,UU);
            shading interp;
            title('GP u');
            
            drawnow;
        end
        %pause;
    end
    toc
    
    %%
    % Final plot
    %
    subplot(2,2,1);
    pcolor(xgrid,tgrid,ks_u_mu);
    shading interp;
    title('KS u');
    colorbar;
    
    subplot(2,2,2);
    pcolor(xgrid,tgrid,UU - ks_u_mu);
    shading interp;
    title('KS u err');
    colorbar;
    
    subplot(2,2,3);
    pcolor(xgrid,tgrid,reshape(mean_pred_u,size(xgrid)));
    shading interp;
    title('GP u');
    colorbar;
    
    subplot(2,2,4);
    pcolor(xgrid,tgrid,UU - reshape(mean_pred_u,size(xgrid)));
    shading interp;
    title('GP u err');
    colorbar;
        
    mse_error = mean((mean_pred_u - u).^2)
    mse_error_ks = mean((ks_u_mu(:) - u).^2)

    %%
    % Final plot 2
    %
    subplot(2,2,1);
    pcolor(xgrid,tgrid,ks_u_mu);
    shading interp;
    title('KS u');
    colorbar;
    
    subplot(2,2,2);
    pcolor(xgrid,tgrid,UU);
    shading interp;
    title('Real u');
    colorbar;
    
    subplot(2,2,3);
    pcolor(xgrid,tgrid,ks_x_mu);
    shading interp;
    title('KS f');
    colorbar;
    
    subplot(2,2,4);
    pcolor(xgrid,tgrid,FF);
    shading interp;
    title('Real f');
    colorbar;
    
    
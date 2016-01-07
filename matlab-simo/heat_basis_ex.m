%
% Solve inverse heat equation in 2d using basis
% function approach.
%

    %%
    % Simulate the data (takes a while)
    %
    heat_sim;

    %%
    % Just use a subset of basis functions 
    %
%    n = 50
    n = size(fun,3);
    fun2 = fun(:,:,1:n);
    lap_e2 = lap_e(1:n);

    %%
    % Spatial covariance function approximation for input
    %
    s = 100;
    ell = 0.05;
    
    se_cov = @(norm_x) s^2 * exp(-norm_x.^2/2/ell^2);
    se_spec = @(norm_w) s^2 * (2*pi) * ell^2 * exp(-ell^2 * norm_w.^2/2);

    Spe = se_spec(sqrt(lap_e2));

    clf;
    subplot(1,2,1);
    pcolor(xx1,xx2,se_cov(sqrt(xx1.^2 + xx2.^2)));
    shading flat
    colorbar;

    N = size(fun2,3);
    
    % Plot k(0,0;x1,x2) = sum_i S(sqrt(L_i)) phi_i(0,0) phi_i(x1,x2)
    fun00 = eigenfun(NN(1:n,:),[0 0]);

    app = zeros(size(xx1));
    for i=1:N
        app = app + Spe(i) * fun00(i) * fun(:,:,i);
    end

    se_spec(sqrt(eigenval(3)))
    Qc = zeros(N);
    for i=1:N
        Qc(i,i) = se_spec(sqrt(eigenval(i)));
    end
    
    subplot(1,2,2);
    pcolor(xx1,xx2,app);
    shading flat
    colorbar;

    %%
    % Temporal part in state-space form
    %
    s_t = 1
    ell_t = 0.5
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
    % Form the state-space model
    %
    
    % evaluator for measurement points
    E = zeros(size(meas_ind,1),N);
    for i=1:size(meas_ind,1)
        for j=1:N
            E(i,j) = fun2(meas_ind(i,1),meas_ind(i,2),j);
        end
    end
    
    % State-space model for the force
    Fgp = kron(eye(N),Fte);
    Lgp = sqrt(qte) * kron(eye(N),Lte);
    Hgpc = kron(eye(N),Hte); % Gives plain coefficients
    Hgp  = kron(E,Hte);      % Projects into measurement points
%    Hgpp = kron(Ep,Hte);     % Projects into prediction points

    Pgp = lyap(Fgp,Lgp*Qc*Lgp');

    % State-space model for the PDE
    D = diffusion;
    lam = decay;
    
    Lpd = eye(N);
    
    Fpd = diag(-decay - diffusion * lap_e2);

    Hpdc = eye(N); % Gives plain coefficients
    Hpd  = E;       % Projects into measurement points
%    Hpdp = Ep;      % Projects into measurement points
 
    % Form the joint state-space model
    F = blkdiag(Fpd,Fgp);
    d = size(Fpd,1);
    F(1:d,d+1:end) = Lpd * Hgpc;
    L = [zeros(d,size(Lgp,2)); Lgp];
    Lc = [Lpd; zeros(size(Lgp,1),size(Lpd,2))];
    H  = [Hpd zeros(size(Hpd,1),size(Fgp,1))];
    Hc = [Hpdc zeros(size(Hpdc,1),size(Fgp,1))];
    Hu  = [zeros(size(Hgp,1),d) Hgp];
    Huc = [zeros(size(Hgpc,1),d) Hgpc];
%    Hp = [Hpdp zeros(size(Hpdp,1),size(Fgp,1))];
%    Hup = [zeros(size(Hgpp,1),d) Hgpp];
    
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
    A = expm(F*dt);
    P = lyap(F,L*Qc*L');
    Q = P - A*P*A';

    fprintf('Done.\n');

    clf;
    subplot(2,2,1);
    spy(A);
    
    subplot(2,2,2);
    spy(Q);

    %%
    % Kalman filter
    %
    clf;
    kf_x_mu = zeros(size(xx1,1),size(xx1,2),length(tt));
    kf_u_mu = zeros(size(xx1,1),size(xx1,2),length(tt));
    
    m = zeros(size(F,1),1);
    %m(1) = x0;
    %m(2) = dx0;
    P = blkdiag(1*eye(size(Fpd,1)),Pgp);
    P0 = P;

    MM = zeros(size(m,1),length(tt));
    PP = zeros(size(P,1),size(P,2),length(tt));
    
    plot_list = 1:10:length(tt);

    R = sd^2*eye(size(YYT,1));
    
    fprintf('Running Kalman filter...\n');
    tic;
    
    min_fc = min(FFT_p(:));
    max_fc = max(FFT_p(:));

    min_uc = min(UUT_p(:));
    max_uc = max(UUT_p(:));
    
    for k=1:length(tt)
        
        m = A*m;
        P = A*P*A' + Q;

        y = YYT(:,k);
        S = H*P*H' + R;
        K = P * H' / S;
        m = m + K * (y - H*m);
        P = P - K * S * K';
        
        MM(:,k) = m;
        PP(:,:,k) = P;
        
        m_c = Hc*m;
        mu_c = Huc*m;
        
        for i=1:size(m_c,1)
            kf_x_mu(:,:,k) = kf_x_mu(:,:,k) + m_c(i) * fun2(:,:,i);
        end
        
        for i=1:size(mu_c,1)
            kf_u_mu(:,:,k) = kf_u_mu(:,:,k) + mu_c(i) * fun2(:,:,i);
        end
        
        if any(plot_list == k)
            subplot(2,2,1);
            pcolor(xx1,xx2,FFT_p(:,:,k));
            shading interp;
            colorbar;
            caxis([min_fc max_fc]);
            title('Real f');
            
            subplot(2,2,3);
            pcolor(xx1,xx2,kf_x_mu(:,:,k));
            shading interp;
            colorbar;
            caxis([min_fc max_fc]);
            title('KF f');
            hold on;
            plot3(meas_x1,meas_x2,YYT(:,k),'x');
            hold off;
        
            subplot(2,2,2);
            pcolor(xx1,xx2,UUT_p(:,:,k));
            shading interp;
            colorbar;
            caxis([min_uc max_uc]);
            title('True u');

            subplot(2,2,4);
            pcolor(xx1,xx2,kf_u_mu(:,:,k));
            shading interp;
            colorbar;
            caxis([min_uc max_uc]);
            title('KF u');
            
            drawnow;
            %pause(0.01);
            %pause;
        end
    end
    toc

    %%
    % Smoother
    %    
    MS = MM;
    PS = PP;
    ms = MM(:,end);
    Ps = PP(:,:,end);
    
    ks_x_mu = zeros(size(kf_x_mu));
    ks_u_mu = zeros(size(kf_u_mu));

    ks_x_mu(:,:,end) = kf_x_mu(:,:,end);
    ks_u_mu(:,:,end) = kf_u_mu(:,:,end);
    
    plot_list = 1:10:length(tt);
    
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
        
        m_c = Hc*ms;
        mu_c = Huc*ms;
        
        for i=1:size(m_c,1)
            ks_x_mu(:,:,k) = ks_x_mu(:,:,k) + m_c(i) * fun2(:,:,i);
        end
        
        for i=1:size(mu_c,1)
            ks_u_mu(:,:,k) = ks_u_mu(:,:,k) + mu_c(i) * fun2(:,:,i);
        end
                
        if any(plot_list == k)
            subplot(2,2,1);
            pcolor(xx1,xx2,FFT_p(:,:,k));
            shading interp;
            colorbar;
            caxis([min_fc max_fc]);
            title('Real f');
            
            subplot(2,2,3);
            pcolor(xx1,xx2,ks_x_mu(:,:,k));
            shading interp;
            colorbar;
            caxis([min_fc max_fc]);
            title('KF f');
            hold on;
            plot3(meas_x1,meas_x2,YYT(:,k),'x');
            hold off;
        
            subplot(2,2,2);
            pcolor(xx1,xx2,UUT_p(:,:,k));
            shading interp;
            colorbar;
            caxis([min_uc max_uc]);
            title('True u');

            subplot(2,2,4);
            pcolor(xx1,xx2,ks_u_mu(:,:,k));
            shading interp;
            colorbar;
            caxis([min_uc max_uc]);
            title('KF u');
            
            drawnow;
            %pause(0.01);
            %pause;
        end
    end
    toc
    
    %%
    % Plot u estimates in the mid times
    %
    k = round(length(tt)/6);
    
    subplot(2,2,1);
    pcolor(xx1,xx2,kf_u_mu(:,:,k));
    shading interp;
    colorbar;
    caxis([min_uc max_uc]);
    title('KF u');

    subplot(2,2,2);
    pcolor(xx1,xx2,ks_u_mu(:,:,k));
    shading interp;
    colorbar;
    caxis([min_uc max_uc]);
    title('KS u');

    subplot(2,2,3);
    pcolor(xx1,xx2,kf_u_mu(:,:,k)-UUT(:,:,k));
    shading interp;
    colorbar;
%    caxis([min_uc max_uc]);
    title('KF u error');

    subplot(2,2,4);
    pcolor(xx1,xx2,ks_u_mu(:,:,k)-UUT(:,:,k));
    shading interp;
    colorbar;
%    caxis([min_uc max_uc]);
    title('KS u error');
    
    %%
    % Plot the maxima (or means) controls at each time
    %
    maxim_r = zeros(1,size(ks_u_mu,3));
    maxim_e = zeros(1,size(ks_u_mu,3));
    for k=1:1:size(ks_u_mu,3)
%        maxim_r(k) = max(max(UUT(:,:,k)));
%        maxim_e(k) = max(max(ks_u_mu(:,:,k)));
        maxim_r(k) = mean(mean(UUT(:,:,k)));
        maxim_e(k) = mean(mean(ks_u_mu(:,:,k)));
    end
    
    clf;
    plot(tt,maxim_r,tt,maxim_e)
    
    %%
    % Plot the maxima (or means) states at each time
    %
    maxim_r = zeros(1,size(ks_x_mu,3));
    maxim_e = zeros(1,size(ks_x_mu,3));
    for k=1:1:size(ks_u_mu,3)
%        maxim_r(k) = max(max(FFT_p(:,:,k)));
%        maxim_e(k) = max(max(ks_x_mu(:,:,k)));
        maxim_r(k) = mean(mean(FFT_p(:,:,k)));
        maxim_e(k) = mean(mean(ks_x_mu(:,:,k)));
    end
    
    clf;
    plot(tt,maxim_r,tt,maxim_e)
    
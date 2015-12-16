%
% Solve Poisson 2D model with state space model and
% basis function expansion.
%

    %%
    % Simulate the data (takes a while)
    %
    %poisson_sim;

    %%
    % Create the state-space model approximation for RBF kernel
    %
    s_t = 1
    ell_t = 0.2

    se_cov_t = @(t) s_t^2 * exp(-t.^2/2/ell_t^2)

%    [Bte,Ate] = se_taylor(4,s_t,ell_t)
    [Bte,Ate] = se_pade(4,8,s_t,ell_t)
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
    s_x = 20
    ell_x = 0.2

    se_cov_x = @(x,y) s_x^2 * exp(-(x-y).^2/2/ell_x^2)
    se_spec_x = @(w) s_x^2 * sqrt(2*pi) * ell_x * exp(-ell_x^2 * w.^2/2);

    N = 50; % How many series terms
    [eigenfun,eigenval,NN] = domain_cartesian(N,1,LL(2));
    lap_e = eigenval(NN);

    tlist = (-LL(1)):dd(1):LL(1);
    xlist = (-LL(2)):dd(2):LL(2);
    [tt,xx] = meshgrid(tlist,xlist);


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
    Spe = se_spec_x(sqrt(lap_e));
    Qc = diag(Spe);

    % Plot k(0,0;x) = sum_i S(sqrt(L_i)) phi_i(0,0) phi_i(x1,x2)
    fun00 = eigenfun(NN,0);

    app = zeros(size(xlist'));
    for i=1:size(NN,1)
        app = app + Spe(i) * fun00(i) * eigenfun(NN(i,:),xlist')
    end

    clf;
    plot(xlist,se_cov_x(xlist,0),xlist,app,'--');
    title('Accuracy of the series approximation of SE');
    legend('Exact SE','Series approx');
    grid on;    

    %%
    % Precompute eigenfunctions
    %
    fun2 = zeros(length(xlist),size(NN,1));

    for k=1:size(NN,1)
        nn = NN(k,:);
        for i=1:length(xlist)
            fun2(i,k) = eigenfun(nn,xlist(i));
        end
    end
    
    clf;
    pcolor(fun2);
    shading flat;
    
    %%
    % Form the joint state-space model
    %

    % Evaluate the basis at measurements
    ut = unique(meas_x(:,1));
    ind = find(meas_x(:,1) == ut(1));
    m_x = meas_x(ind,2);
    
    E = zeros(length(m_x),N);
    for i=1:length(m_x)
        for j=1:N
            tx = m_x(i);
            E(i,j) = eigenfun(j,tx);
        end
    end

    % Evaluate the basis at prediction points
    Ep = zeros(length(xlist),N);
    for i=1:length(xlist)
        for j=1:N
            tx = xlist(i);
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
    % dvf/dt = [0 1; nabla^2  -2 sqrt(-nabla^2)] vf + [0;-1] u
    Fpd = zeros(2*N);
    Lpd = kron(eye(N),[0;1]);
    
    for i=1:N
        Fpd(2*i-1,2*i) = 1;
        Fpd(2*i,2*i-1) = -lap_e(i);
        Fpd(2*i,2*i)   = -2*sqrt(lap_e(i));
    end

    Hpdc = kron(eye(N),[1 0]); % Gives plain coefficients
    Hpd  = kron(E,[1 0]);      % Projects into measurement points
    Hpdp = kron(Ep,[1 0]);     % Projects into prediction points
     
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
    dt = dd(1);
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
    % Check that the forward model works
    %
    uu_c2 = zeros(size(NN,1),length(tlist));
    uu_p2 = zeros(size(uu));
    for k=1:length(tlist)
        for i=1:size(NN,1)
            uu_c2(i,k) = myinner(uu(k,:)',fun2(:,i));
            uu_p2(k,:) = uu_p(k,:) + uu_c2(i,k) * fun2(:,i)';
        end
    end

    Apd = expm(Fpd*dt);
    %Bpd = Lpd*dt;
    tmp = expm([Fpd Lpd; zeros(size(Lpd')) zeros(size(Lpd,2))]*dt);
    n = size(Fpd,1);
    Bpd = tmp(1:n,n+1:end);
    
    clf;
    m = zeros(size(Fpd,1),1);
    ff_fwd = zeros(size(uu));
    for k=1:size(uu,1)
        m = Apd*m + Bpd*uu_c2(:,k);
        ff_fwd(k,:) = (Hpdp*m)';

        subplot(1,2,1);
        pcolor(ff_p)
        shading flat;
        colorbar;
        ac = caxis;

        subplot(1,2,2);
        pcolor(ff_fwd)
        shading flat;
        colorbar;
%        caxis(ac);
        
        drawnow;

%        pause(0.01);
    
    end
    
    %%
    % Kalman filter
    %
    clf;
    m_t = unique(meas_x(:,1));
    ind = find(meas_x(:,1) == ut(1));
    m_x = meas_x(ind,2);

    kf_x_mu = zeros(length(xlist),length(tlist));
    kf_u_mu = zeros(length(xlist),length(tlist));
    
    m = zeros(size(F,1),1);
    P = blkdiag(0*eye(size(Fpd,1)),Pgp);
    P0 = P;

    MM = zeros(size(m,1),length(tlist));
    PP = zeros(size(P,1),size(P,2),length(tlist));
    
    fprintf('Running Kalman filter...\n');

    tic;
    plot_list = 1:10:length(tlist);

    for k=1:length(tlist)
        
        m = A*m;
        P = A*P*A' + Q;

        MM(:,k) = m;
        PP(:,:,k) = P;

        kf_x_mu(:,k) = Hp*m;
        kf_u_mu(:,k) = Hup*m;

        ind = find(abs(meas_x(:,1) - tlist(k)) < eps);
        if ~isempty(ind)
            y = Y(ind)
        
            R = sd^2 * eye(length(y));
            S = H*P*H' + R;
            K = P * H' / S;
            m = m + K * (y - H*m);
            P = P - K * S * K';
        
            MM(:,k) = m;
            PP(:,:,k) = P;
        end
        
        kf_x_mu(:,k) = Hp*m;
        kf_u_mu(:,k) = Hup*m;

        if any(plot_list == k)
            subplot(2,2,1);
            pcolor(xx,tt,kf_x_mu);
            shading interp;
            title('KF f');
            
            subplot(2,2,2);
            pcolor(xx,tt,ff_p);
            shading interp;
%            title('GP f');
            title('Actual f');
            
            subplot(2,2,3);
            pcolor(xx,tt,kf_u_mu);
            shading interp;
            title('KF u');
            
            subplot(2,2,4);
            pcolor(xx,tt,uu);
            shading interp;
%            title('GP u');
            title('Actual u');

            drawnow;
        end

        %pause;
    end
 
    % Dirichlet boundary condition at end
    R = 1e-6 * eye(size(Hp,1));
    S = Hp*P*Hp' + R;
    K = P * Hp' / S;
    m = m + K * (-Hp*m);
    P = P - K * S * K';

    MM(:,k) = m;
    PP(:,:,k) = P;
        
    kf_x_mu(:,k) = Hp*m;
    kf_u_mu(:,k) = Hup*m;
        
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
    
    plot_list = 1:10:size(MM,2);
    
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
            pcolor(xx,tt,ks_x_mu);
            shading interp;
            title('KS f');
            
            subplot(2,2,2);
            pcolor(xx,tt,ff_p);
            shading interp;
%            title('GP f');
            title('Actual f');
            
            subplot(2,2,3);
            pcolor(xx,tt,ks_u_mu);
            shading interp;
            title('KS u');
            
            subplot(2,2,4);
            pcolor(xx,tt,uu);
            shading interp;
%            title('GP u');
            title('Actual u');
            
            drawnow;
        end
        %pause;
    end
    
    toc
    
    
%
% Test KF/RTS with a SE latent force model
%

    %%
    % Simulate
    %
    quiet = 1;
    clf;
    anim_open;

    %%
    % Optimize
    %
    opt = optimset('Display','iter');
    s_ell0 = [1 5] + rand;
    s_ell = fminsearch(@(s_ell) spring_se_optfun(s_ell,sd,b,order,dt,Y),s_ell0,opt)
    
    %%
    % Create the state-space model approximation for RBF kernel
    %
    addpath ..
%    s = 1; % Scale
%    ell = 5;
    s = s_ell(1);
    ell = s_ell(2);
    
    order = 4; % Order of state-space model
    
    se_cov  = @(t) s^2 * exp(-t.^2/2/ell^2);

%    [Bgp,Agp] = se_taylor(order,s,ell);
    [Bgp,Agp] = se_pade(4,8,s,ell)
    [Fgp,Lgp,qgp,Hgp] = ratspec_to_ss(Bgp,Agp);
    [Fgp,Lgp,Hgp] = ss_balance(Fgp,Lgp,Hgp);

    Fgp
    Lgp
    Hgp
    qgp
    
    tau = -5*ell:ell/10:5*ell;
    c = ss_cov(tau,Fgp,Lgp,qgp,Hgp);
    cond(Fgp)

    Pgp = lyapchol(Fgp,Lgp*sqrt(qgp));
    Pgp = Pgp' * Pgp;
    
    clf;
    plot(tau,se_cov(tau),tau,c,'r--');
    title('Accuracy of the state-space approximation of SE');
    legend('Exact SE','State-space Approx');
    grid on;


    %%
    % Form the joint model
    %
    Fsp = [0 1; -1 -b];
    Lsp = [0;1];
    Hsp = [1 0];
    
    Fjm = blkdiag(Fsp,Fgp);
    d = size(Fsp,1);
    Fjm(1:d,d+1:end) = Lsp * Hgp
    Ljm  = [zeros(d,1); Lgp]
    Ljmc = [Lsp; zeros(size(Lgp,1),1)];
    Hjm  = [Hsp zeros(1,size(Fgp,1))]
    Hjmu = [zeros(1,d) Hgp]
    qjm = qgp;
    R = sd^2;
    
    cond(Fjm)
    
    %%
    % Run Kalman filter
    %
    kf_x_mu = [];
    kf_x_V  = [];

    kf_u_mu = [];
    kf_u_V  = [];
    
    m = zeros(size(Fjm,1),1);
    P = blkdiag(eye(size(Fsp,1)),Pgp)
    
    MM = zeros(size(m,1),length(T));
    PP = zeros(size(P,1),size(P,2),length(T));
    
    [Ajm,Qjm] = lti_disc(Fjm,Ljm,qjm,dt);
    
    nlog_lh = 0;
    
    for k=1:length(Y)
        m = Ajm*m;
        P = Ajm*P*Ajm' + Qjm;
        
        if ~isnan(Y(k))
            S = Hjm*P*Hjm' + R;
            K = P * Hjm' / S;
            v = (Y(k) - Hjm*m);
            m = m + K * v;
            P = P - K * S * K';

            nlog_lh = nlog_lh + 0.5*log(S) + 0.5*v^2/S;
        end
        
        MM(:,k) = m;
        PP(:,:,k) = P;
        
        kf_x_mu = [kf_x_mu Hjm*m];
        kf_x_V  = [kf_x_V  Hjm*P*Hjm'];

        kf_u_mu = [kf_u_mu Hjmu*m];
        kf_u_V  = [kf_u_V  Hjmu*P*Hjmu'];
        
    end
    
    ind = find(~isnan(Y));
    
    subplot(2,1,1);
    plot(T,kf_x_mu,T(ind),Y(ind),'.')

    subplot(2,1,2);
    plot(T,kf_u_mu,'--',T,u_ext)
    
    %%
    % Smoother
    %    
    MS = MM;
    PS = PP;
    ms = MM(:,end);
    Ps = PP(:,:,end);
    
    ks_x_mu = Hjm*ms;
    ks_x_V  = Hjm*Ps*Hjm';

    ks_u_mu = Hjmu*ms;
    ks_u_V  = Hjmu*Ps*Hjmu';
    
    for k=size(MM,2)-1:-1:1
        m = MM(:,k);
        P = PP(:,:,k);
        
        mp = Ajm*m;
        Pp = Ajm*P*Ajm' + Qjm;
        
        G = P*Ajm'/Pp;
        ms = m + G * (ms - mp);
        Ps = P + G * (Ps - Pp) * G';
        
        MS(:,k) = ms;
        PS(:,:,k) = Ps;
        
        ks_x_mu = [Hjm*ms ks_x_mu];
        ks_x_V  = [Hjm*Ps*Hjm' ks_x_V];

        ks_u_mu = [Hjmu*ms ks_u_mu];
        ks_u_V  = [Hjmu*Ps*Hjmu' ks_u_V];        
    end
    
    ind = find(~isnan(Y));

    subplot(2,1,1);
    plot(T,ks_x_mu,T(ind),Y(ind),'.')

    subplot(2,1,2);
    plot(T,ks_u_mu,'r--',T,u_ext)
    
    %%
    % Test Basic LQ control with the full LFM
    %
    G = lqr(Fsp,Lsp,diag([1 1]),0.05);
    
    kf_x_mu1 = [];
    kf_x_V1  = [];

    kf_u_mu1 = [];
    kf_u_V1  = [];
    
    x = [-1;0];
    Z1 = [];
    
    m = zeros(size(Fjm,1),1);
    P = blkdiag(eye(size(Fsp,1)),Pgp)
    
    MM1 = zeros(size(m,1),length(T));
    PP1 = zeros(size(P,1),size(P,2),length(T));
    
    [Ajm,Qjm] = lti_disc(Fjm,Ljm,qjm,dt);
    
    rng(1);

    for k=1:length(Y)
        u = - G*m(1:2);
        
        x = x + dt*(Fsp*x + Lsp*u + Lsp*u_ext(k));        
        
        if rem(k,meas_step) == 0
            z = x(1) + sd*randn;
        else
            z = NaN;
        end
        
        Z1 = [Z1 z];
        
        m = Ajm*m + Ljmc*u*dt;
        P = Ajm*P*Ajm' + Qjm;
        
        if ~isnan(z)
            S = Hjm*P*Hjm' + R;
            K = P * Hjm' / S;
            m = m + K * (z - Hjm*m);
            P = P - K * S * K';
        end
        
        MM1(:,k) = m;
        PP1(:,:,k) = P;
        
        kf_x_mu1 = [kf_x_mu1 Hjm*m];
        kf_x_V1  = [kf_x_V1  Hjm*P*Hjm'];

        kf_u_mu1 = [kf_u_mu1 Hjmu*m];
        kf_u_V1  = [kf_u_V1  Hjmu*P*Hjmu'];
        
    end
    
    ind = find(~isnan(Z1));
    
    subplot(2,1,1);
    plot(T,kf_x_mu1,T(ind),Z1(ind),'.')
    grid on;

    subplot(2,1,2);
    plot(T,kf_u_mu1,'--',T,u_ext)
    grid on;
    
    
    %%
    % Test LQ control with the full LFM
    %
    G = lqr(Fjm,Ljmc,eye(size(Fjm,1)),0.05);
    
    kf_x_mu2 = [];
    kf_x_V2  = [];

    kf_u_mu2 = [];
    kf_u_V2  = [];
    
    x = [-1;0];
    Z2 = [];
    
    m = zeros(size(Fjm,1),1);
    P = blkdiag(eye(size(Fsp,1)),Pgp)
    
    MM2 = zeros(size(m,1),length(T));
    PP2 = zeros(size(P,1),size(P,2),length(T));
    
    [Ajm,Qjm] = lti_disc(Fjm,Ljm,qjm,dt);
    
    rng(1);
    
    for k=1:length(Y)
        u = - G*m;
        
        x = x + dt*(Fsp*x + Lsp*u + Lsp*u_ext(k));        
        
        if rem(k,meas_step) == 0
            z = x(1) + sd*randn;
        else
            z = NaN;
        end
        
        Z2 = [Z2 z];
        
        m = Ajm*m + Ljmc*u*dt;
        P = Ajm*P*Ajm' + Qjm;
        
        if ~isnan(z)
            S = Hjm*P*Hjm' + R;
            K = P * Hjm' / S;
            m = m + K * (z - Hjm*m);
            P = P - K * S * K';
        end
        
        MM2(:,k) = m;
        PP2(:,:,k) = P;
        
        kf_x_mu2 = [kf_x_mu2 Hjm*m];
        kf_x_V2  = [kf_x_V2  Hjm*P*Hjm'];

        kf_u_mu2 = [kf_u_mu2 Hjmu*m];
        kf_u_V2  = [kf_u_V2  Hjmu*P*Hjmu'];
        
    end
    
    ind = find(~isnan(Z2));
    
    subplot(2,1,1);
    
    plot(T,kf_x_mu2,T(ind),Z2(ind),'.')
    grid on;

    subplot(2,1,2);
    plot(T,kf_u_mu2,'--',T,u_ext)
    grid on;
    
    
    
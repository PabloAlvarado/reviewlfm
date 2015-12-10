%
% State space RBF LFM ODE demo
%

    %%
    % Load the data
    %
    load ../matlab-mauricio/dataODEwithRBF.mat;
    
    %%
    % Create the state-space model approximation for RBF kernel
    %
    s = 1; % Scale
    ell = lengthscale;
    
    order = 12; % Order of state-space model
    
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
    % State space representation of the spring model:
    %
    %   x'' + lam x' + gam x = sen u
    %
    %  [x;x']' = [0 1; -gam -lam] + [0;sen] u
    %  y = [1 0] [x;x']
    %
    lam = damperc;
    gam = springc;
    sen = sensitivityc;
    
    Fsp = [0 1; -gam -lam]
    Lsp = [0;sen]
    Hsp = [1 0]

    cond(Fsp)
    
    %%
    % Create the joint state-space
    %
    F = blkdiag(Fsp,Fgp);
    d = size(Fsp,1);
    F(1:d,d+1:end) = Lsp * Hgp
    L = [zeros(d,1); Lgp]
    H = [Hsp zeros(1,size(Fgp,1))]
    Hu = [zeros(1,d) Hgp]
    q = qgp;
    R = sigma2;
    
%    [F,L,H,T] = ss_balance(F,L,H);
%    Hu = Hu * T;
    cond(F)
    
    %%
    % Run Kalman filter and smoother
    %
    kf_x_mu = [];
    kf_x_V  = [];

    kf_u_mu = [];
    kf_u_V  = [];
    
    m = zeros(size(F,1),1);
    %m(1) = x0;
    %m(2) = dx0;
    P = blkdiag(0*eye(size(Fsp,1)),Pgp)
    
    MM = zeros(size(m,1),length(t));
    PP = zeros(size(P,1),size(P,2),length(t));
    
    dt = t(2) - t(1);
    [A,Q] = lti_disc(F,L,q,dt);
    
    count = 0;
    yy = nan * ones(size(f));
    yy(indexk) = yk;
        
    for k=1:length(yy)
        if k > 1
            m = A*m;
            P = A*P*A' + Q;
        end
        
        if ~isnan(yy(k))
            S = H*P*H' + R;
            K = P * H' / S;
            m = m + K * (yy(k) - H*m);
            P = P - K * S * K';
            
        end
        
        MM(:,k) = m;
        PP(:,:,k) = P;
        
        kf_x_mu = [kf_x_mu H*m];
        kf_x_V  = [kf_x_V  H*P*H'];

        kf_u_mu = [kf_u_mu Hu*m];
        kf_u_V  = [kf_u_V  Hu*P*Hu'];
        
    end
    
    subplot(2,1,1);
    plot(t,kf_x_mu,t,f,tk,yk,'o')

    subplot(2,1,2);
    plot(t,kf_u_mu,t,mean_pred_u,'--')

    %%
    % Smoother
    %    
    MS = MM;
    PS = PP;
    ms = MM(:,end);
    Ps = PP(:,:,end);
    
    ks_x_mu = H*ms;
    ks_x_V  = H*Ps*H';

    ks_u_mu = Hu*ms;
    ks_u_V  = Hu*Ps*Hu';
    
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
        
        ks_x_mu = [H*ms ks_x_mu];
        ks_x_V  = [H*Ps*H' ks_x_V];

        ks_u_mu = [Hu*ms ks_u_mu];
        ks_u_V  = [Hu*Ps*Hu' ks_u_V];        
    end
    
    subplot(3,1,1);
    plot(t,ks_x_mu,t,f,tk,yk,'o')

    subplot(3,1,2);
    plot(t,ks_u_mu,t,mean_pred_u,'r--')
    
    subplot(3,1,3);
    %semilogy(t,ks_u_V,t,var_pred_u,'r--')
    plot(t,ks_u_V,t,var_pred_u,'r--')
    
    %%
    % An ODE test
    %
    sol = zeros(size(f));
    dt = t(2) - t(1);    
    xx = [0;0];

    for k=1:length(sol)
        sol(k) = xx(1);
        xx = xx + (Fsp*xx) * dt + Lsp * u(k) * dt;
    end
    %plot(t,sol,t,f);
        
    %%
    % Plot the results
    %

    mu = ks_u_mu;
    V  = ks_u_V;
    
    clf;
    p=patch([t' fliplr(t')], ...
            [mu + 1.96*sqrt(V) fliplr(mu - 1.96*sqrt(V))],1);
    set(p,'EdgeColor','none')
    colormap(0.8*[1 1 1]);
    hold on;
    h = plot(t,u,'k--',t,mu,'b-',tk,yk,'ro');
    set(h,'Linewidth',3);
    set(h(1),'Linewidth',2);
    set(h,'MarkerSize',10);
    grid on;
%    axis([0 11 -2.5 2]);
    legend(h,'True function','Regression mean','Observations',3);

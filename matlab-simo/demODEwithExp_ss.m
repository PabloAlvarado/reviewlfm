%
% State space Matern LFM ODE demo
%

    %%
    % Reset everything
    %
    clc
    clear
    close all
    s = RandStream('mt19937ar', 'Seed', 1);
    RandStream.setDefaultStream(s);
%    RandStream.setGlobalStream(s);
    
    damperc = 40; % Damper constant = lambda if mass = 1
    springc = 10; % Spring constant = gamma  if mass = 1
    sensitivityc = 100; % Sensitivity parameter


    %%
    % Exp in state-space form
    %
    magnSigma2 = 1;
    lengthScale = 0.1;

    matern = @(tau) magnSigma2 * exp(-abs(tau)./lengthScale);
    
    % Derived constants
    lambda = 1/lengthScale;
  
    % Feedback matrix
    Fgp = -lambda

    % Noise effect matrix
    Lgp = 1

    % Spectral density
    qgp = 2*lambda*magnSigma2

    % Observation model
    Hgp = 1
  
    tau = -5*lengthScale:lengthScale/10:5*lengthScale;
    c = ss_cov(tau,Fgp,Lgp,qgp,Hgp);
    cond(Fgp)

    Pgp = lyapchol(Fgp,Lgp*sqrt(qgp));
    Pgp = Pgp' * Pgp;
    
    clf;
    plot(tau,matern(tau),tau,c,'r--');
    title('Accuracy of the state-space approximation of SE');
    legend('Exact Exp','State-space Approx');
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
    
%    [F,L,H,T] = ss_balance(F,L,H);
%    Hu = Hu * T;
    cond(F)
        
    
    %%
    % Generate data
    %
    N = 100; % Number of time points
    t = linspace(0, 1, N)'; 

    f = zeros(size(t));
    u = zeros(size(t));

    m = zeros(size(F,1),1);
    %m(1) = x0;
    %m(2) = dx0;
    P = blkdiag(0*eye(size(Fsp,1)),Pgp);
    C = blkdiag(0*eye(size(Fsp,1)),chol(Pgp,'lower'));
    
    x = m + C*randn(size(F,1),1);

    dt = t(2) - t(1);
    [A,Q] = lti_disc(F,L,q,dt);

    C = chol(Q,'lower');
    
    for k=1:length(t)
        if k > 1
            x = A * x + C * randn(size(F,1),1);
        end
        f(k) = H * x;
        u(k) = Hu * x;
    end
    
    Nk = 30;
    index = randperm(N);
    indexk = sort(index(1:Nk));
    fk = f(indexk);
    tk = t(indexk);
    % We add some noise to the measurements fk,
    sigma2 = 0.01*var(f);
    yk = fk +  sqrt(sigma2)*randn(Nk, 1);

    %%
    % Run Kalman filter and smoother
    %
    R = sigma2;
    
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
    %plot(t,kf_u_mu,t,mean_pred_u,'--')
    plot(t,kf_u_mu,t,u,'--')

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
    %plot(t,ks_u_mu,t,mean_pred_u,'r--')
    plot(t,ks_u_mu,t,u,'r--')
    
    subplot(3,1,3);
    %semilogy(t,ks_u_V,t,var_pred_u,'r--')
    plot(t,ks_u_V)
    
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

    clf;
    subplot(2,1,1);
    mu = ks_x_mu;
    V  = ks_x_V;
    
    p=patch([t' fliplr(t')], ...
            [mu + 1.96*sqrt(V) fliplr(mu - 1.96*sqrt(V))],1);
    set(p,'EdgeColor','none')
    colormap(0.8*[1 1 1]);
    hold on;
    h = plot(t,f,'k--',t,mu,'b-',tk,yk,'ro');
    set(h,'Linewidth',3);
    set(h(1),'Linewidth',2);
    set(h,'MarkerSize',10);
    grid on;
    legend(h,'True function','Regression mean','Observations',3);

    subplot(2,1,2);
    mu = ks_u_mu;
    V  = ks_u_V;
    
    p=patch([t' fliplr(t')], ...
            [mu + 1.96*sqrt(V) fliplr(mu - 1.96*sqrt(V))],1);
    set(p,'EdgeColor','none')
    colormap(0.8*[1 1 1]);
    hold on;
    h = plot(t,u,'k--',t,mu,'b-');
    set(h,'Linewidth',3);
    set(h(1),'Linewidth',2);
    set(h,'MarkerSize',10);
    grid on;
    legend(h,'True force','Regression mean',3);
    
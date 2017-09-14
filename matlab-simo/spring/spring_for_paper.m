%
% Generate spring results for the paper 
%

    %%
    % Random seed and parameters
    %
    s = RandStream('mt19937ar', 'Seed', 1e8);
    RandStream.setGlobalStream(s);
    
    spring_param;

    s_ell = [1.4173; 2.9544];
    
    %%
    % Create the state-space model approximation for RBF kernel
    %
%    s = 1; % Scale
%    ell = 5;
    s = s_ell(1);
    ell = s_ell(2);
    
    se_cov  = @(t) s^2 * exp(-t.^2/2/ell^2);

%    order = 6;
%    [Bgp,Agp] = se_taylor(order,s,ell);
    [Bgp,Agp] = se_pade(4,8,s,ell)
    [Fgp,Lgp,qgp,Hgp] = ratspec_to_ss(Bgp,Agp);
    [Fgp,Lgp,Hgp] = ss_balance(Fgp,Lgp,Hgp);

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
    % Run basic LQ
    %    
    [G,S_lq,E_lq] = lqr(Fsp,Lsp,diag([1 1]),0.05);
    
    T = tt;
    
    sd_b = 0.01;
    R_b = sd_b^2;
    
    m = zeros(size(Fjm,1),1);
    P = blkdiag(eye(size(Fsp,1)),Pgp)
    
    MM1 = zeros(size(m,1),length(T));
    PP1 = zeros(size(P,1),size(P,2),length(T));
    UU1 = zeros(1,length(T));
    
    [Ajm,Qjm] = lti_disc(Fjm,Ljm,qjm,dt);
    
    n_mc = 1;

    for iter=1:n_mc    
        x = [-1;0];
    
        kf_x_mu1 = [];
        kf_x_V1  = [];

        kf_u_mu1 = [];
        kf_u_V1  = [];
    
        Z1 = [];
        XR1 = [];
        
        for k=1:length(T)
            u = - G*m(1:2);
            UU1(k) = u;

            x = x + dt*(Fsp*x + Lsp*u + Lsp*u_ext(k));        

            if rem(k,meas_step) == 0
                z = x(1) + sd_b*randn;
            else
                z = NaN;
            end

            Z1 = [Z1 z];

            m = Ajm*m + Ljmc*u*dt;
            P = Ajm*P*Ajm' + Qjm;

            if ~isnan(z)
                S = Hjm*P*Hjm' + R_b;
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

            XR1 = [XR1 x];
        end
        
        err1 = sqrt(mean(XR1(1,:).^2))
%        pause;
    end
    
    ind = find(~isnan(Z1));
    
    clf;
    plot(T,XR1(1,:))
    grid on;
    title('Position');

    %%
    % LFM LQ control
    %
    [G,S_lfm,E_lfm] = lqr(Fjm,Ljmc,diag([1 1 zeros(1,size(Fjm,1)-2)]),0.05);
    
    m = zeros(size(Fjm,1),1);
    P = blkdiag(eye(size(Fsp,1)),Pgp);
    
    MM2 = zeros(size(m,1),length(T));
    PP2 = zeros(size(P,1),size(P,2),length(T));
    UU2 = zeros(1,length(T));
    
    [Ajm,Qjm] = lti_disc(Fjm,Ljm,qjm,dt);
    
    for iter=1:n_mc    
        x = [-1;0];
    
        kf_x_mu2 = [];
        kf_x_V2  = [];

        kf_u_mu2 = [];
        kf_u_V2  = [];

        Z2 = [];
        
        XR2 = [];

        for k=1:length(T)
            u = - G*m;
            UU2(k) = u;

            x = x + dt*(Fsp*x + Lsp*u + Lsp*u_ext(k));        

            if rem(k,meas_step) == 0
                z = x(1) + sd_b*randn;
            else
                z = NaN;
            end

            Z2 = [Z2 z];

            m = Ajm*m + Ljmc*u*dt;
            P = Ajm*P*Ajm' + Qjm;

            if ~isnan(z)
                S = Hjm*P*Hjm' + R_b;
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
            
            XR2 = [XR2 x];

        end
        
        err2 = sqrt(mean(XR2(1,:).^2))
%        pause;
    end
    
    ind = find(~isnan(Z2));
    
    clf;
    plot(T,XR2(1,:))
    grid on;
    title('Position');

    %%
    % Plot both results
    %
    h = plot(T,XR1(1,:),'--',T,XR2(1,:))
    set(h(1),'LineWidth',2);
    set(h(1),'Color',[0.7 0.7 0.7]);
    set(h(2),'LineWidth',1);
    set(h(2),'Color',[0.0 0.0 0.0]);
    grid on;
    legend('Basic LQR','LFM LQR','Location','SE');
    xlabel('Time{\it t}');
    ylabel('Position{\it f(t)}');
    
    print -depsc cntl_spring.eps;
    
    %%
    %
    %
    clf;
    h = plot(T,kf_u_mu2,'--',T,u_ext)
    set(h(1),'LineWidth',2);
    set(h(1),'Color',[0.7 0.7 0.7]);
    set(h(2),'LineWidth',1);
    set(h(2),'Color',[0.0 0.0 0.0]);
    grid on;
    legend('Kalman estimate','True force','Location','SE');
    xlabel('Time{\it t}');
    ylabel('Force{\it u(t)}');
    
    print -depsc cntl_force.eps;
    
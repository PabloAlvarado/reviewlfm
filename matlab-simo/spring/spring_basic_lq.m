%
% Test Basic LQ control with the full LFM
%
    
    %%
    % Random seed and parameters
    %
    s = RandStream('mt19937ar', 'Seed', 1e8);
    RandStream.setGlobalStream(s);
    
    %%
    % Run the LQ controller
    %
    [G,S_lq,E_lq] = lqr(Fsp,Lsp,diag([1 1]),0.05);
    
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
    UU1 = zeros(1,length(T));
    
    [Ajm,Qjm] = lti_disc(Fjm,Ljm,qjm,dt);
    
    rng(1);

    for k=1:length(Y)
        u = - G*m(1:2);
        UU1(k) = u;
        
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
    title('Position');

    subplot(2,1,2);
    plot(T,kf_u_mu1,'--',T,u_ext)
    grid on;
    title('Latent force');
    
    %%
    % Plot control
    %
    subplot(2,1,1);
    plot(T,UU1);
    grid on;
    title('Control');

    subplot(2,1,2);
    plot(T,kf_u_mu1+UU1,'--',T,u_ext+UU1)
    grid on;
    title('Latent force + control');

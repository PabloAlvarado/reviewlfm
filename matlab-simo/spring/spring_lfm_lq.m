%
% Test LQ control with the full LFM
%
    
    %%
    % Check the theory on the relation to LQ solution
    %
    U = 0.05;
    X = diag([1 1]);
    [G_lq,S_lq,E_lq] = lqr(Fsp,Lsp,X,U);
    [G_lfm,S_lfm,E_lfm] = lqr(Fjm,Ljmc,blkdiag(X,zeros(size(Fjm,1)-2)),U);
    
    K12 = sylvester(S_lq*Lsp/U*Lsp'-Fsp',-Fgp,S_lq*Lsp*Hgp)
    K22 = lyap(-Fgp',K12'*Lsp/U*Lsp'*K12 - Hgp'*Lsp'*K12 - K12'*Lsp*Hgp)
    S_lfm2 = [S_lq K12; K12' K22];
    
    S_norm = norm(S_lfm-S_lfm2)
    G_lfm2 = [U\Lsp'*S_lq U\Lsp'*K12];
    
    G_norm = norm(G_lfm-G_lfm2)
 
    %%
    % LQ control
    %
    [G,S_lfm,E_lfm] = lqr(Fjm,Ljmc,diag([1 1 zeros(1,size(Fjm,1)-2)]),0.05);
    
    kf_x_mu2 = [];
    kf_x_V2  = [];

    kf_u_mu2 = [];
    kf_u_V2  = [];
    
    x = [-1;0];
    Z2 = [];
    
    m = zeros(size(Fjm,1),1);
    P = blkdiag(eye(size(Fsp,1)),Pgp);
    
    MM2 = zeros(size(m,1),length(T));
    PP2 = zeros(size(P,1),size(P,2),length(T));
    UU2 = zeros(1,length(T));
    
    [Ajm,Qjm] = lti_disc(Fjm,Ljm,qjm,dt);
    
    rng(1);
    
    for k=1:length(Y)
        u = - G*m;
        UU2(k) = u;
        
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
    title('Position');

    subplot(2,1,2);
    plot(T,kf_u_mu2,'--',T,u_ext)
    grid on;
    title('Latent force');
    
    %%
    % Plot control
    %
    subplot(2,1,1);
    plot(T,UU2);
    grid on;
    title('Control');

    subplot(2,1,2);
    plot(T,kf_u_mu2+UU2,'--',T,u_ext+UU2)
    grid on;
    title('Latent force + control');
    

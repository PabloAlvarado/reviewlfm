%
% GP Matern LFM ODE demo with covariances calculated
% via a state-space method
%

    %%
    % Reset everything
    %
    load dataODEwithMatern;

    %%
    % Matern in state-space form
    %
    magnSigma2 = 1;
    lengthScale = 0.1;

    matern = @(tau) magnSigma2 * ...
    exp(-sqrt(3)*abs(tau)./lengthScale).*(sqrt(3)*abs(tau)/lengthScale+1);
    
    % Derived constants
    lambda = sqrt(3)/lengthScale;
  
    % Feedback matrix
    Fgp = [0,          1;
           -lambda^2,  -2*lambda]

    % Noise effect matrix
    Lgp = [0; 1]

    % Spectral density
    qgp = 12*sqrt(3)/lengthScale^3*magnSigma2

    % Observation model
    Hgp = [1 0]
  
    tau = -5*lengthScale:lengthScale/10:5*lengthScale;
    c = ss_cov(tau,Fgp,Lgp,qgp,Hgp);
    cond(Fgp)

    Pgp = lyapchol(Fgp,Lgp*sqrt(qgp));
    Pgp = Pgp' * Pgp;
    
    clf;
    plot(tau,matern(tau),tau,c,'r--');
    title('Accuracy of the state-space approximation of SE');
    legend('Exact Matern','State-space Approx');
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
    P0 = blkdiag(0*eye(size(Fsp,1)),Pgp)

    %%
    % Compute the posteriors of u and x (i.e. f)
    %
    Kfy = zeros(length(f),length(yk));
    Kuy = zeros(length(u),length(yk));
    Kyy = zeros(length(yk),length(yk));
    for j=1:length(yk)
        for i=1:length(u)
            Kfy(i,j) = H*ss_crosscov(t(i),tk(j),F,L,q,P0)*H';
            Kuy(i,j) = Hu*ss_crosscov(t(i),tk(j),F,L,q,P0)*H';
        end
        for i=1:length(yk)
            Kyy(i,j) = H*ss_crosscov(tk(i),tk(j),F,L,q,P0)*H';
            if i == j
                Kyy(i,j) = Kyy(i,j) + sigma2;
            end
        end
    end
    Kuu = zeros(length(u),length(u));
    Kff = zeros(length(f),length(f));
    for j=1:length(u)
        for i=1:length(u)
            Kff(i,j) = H*ss_crosscov(t(i),t(j),F,L,q,P0)*H';
            Kuu(i,j) = Hu*ss_crosscov(t(i),t(j),F,L,q,P0)*Hu';
        end
    end
    gp_u_mu = Kuy / Kyy * yk;
    gp_u_V  = diag(Kuu - Kuy / Kyy * Kuy');
    gp_x_mu = Kfy / Kyy * yk;
    gp_x_V  = diag(Kff - Kfy / Kyy * Kfy');
    
    %%
    % Plot the results
    %

    clf;
    subplot(2,1,1);
    mu = gp_x_mu';
    V  = gp_x_V';
    
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
    mu = gp_u_mu';
    V  = gp_u_V';
    
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
    
    %%
    % Compare KF and GP solutions
    %
    clf;
    subplot(2,1,1);
    plot(t,gp_u_mu,'b',t,ks_u_mu,'r--');
    hold on;
    plot(t,gp_x_mu,'k',t,ks_x_mu,'g--');
    
    subplot(2,1,2);
    semilogy(t,gp_u_V,'b',t,ks_u_V,'r--');
    hold on;
    semilogy(t,gp_x_V,'k',t,ks_x_V,'g--');
    
    
    
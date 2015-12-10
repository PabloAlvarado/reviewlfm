%
% Solve Poisson 2D model with basis function expansion.
%

    %%
    % Simulate the data (takes a while)
    %
    poisson_sim;

    %%
    % Just use a subset of basis functions 
    %
    n = 100;
    fun2 = fun(:,:,1:n);
    lap_e2 = lap_e(1:n);


    %%
    % Covariance function approximation for input
    %
    ell = 0.2;
    s = 20;

    se_cov = @(norm_x) s^2 * exp(-norm_x.^2/2/ell^2);
    se_spec = @(norm_w) s^2 * (2*pi) * ell^2 * exp(-ell^2 * norm_w.^2/2);

    Spe = se_spec(sqrt(lap_e2));

    clf;
    subplot(1,2,1);
    pcolor(xx1,xx2,se_cov(sqrt(xx1.^2 + xx2.^2)));
    shading flat
    colorbar;

    % Plot k(0,0;x1,x2) = sum_i S(sqrt(L_i)) phi_i(0,0) phi_i(x1,x2)
    fun00 = eigenfun(NN(1:n,:),[0 0]);

    app = zeros(size(xx1));
    for i=1:size(fun2,3)
        app = app + Spe(i) * fun00(i) * fun(:,:,i);
    end

    subplot(1,2,2);
    pcolor(xx1,xx2,app);
    shading flat
    colorbar;

    %%
    % Covariance function of the output
    %
    app = zeros(size(xx1));
    for i=1:size(fun2,3)
        app = app + Spe(i)/lap_e2(i)^2 * fun00(i) * fun(:,:,i);
    end

    clf;
    pcolor(xx1,xx2,app);
    shading flat
    colorbar;


    %%
    % Solve the input with GP regression
    % y = C*G*u + e
    %
    P = diag(Spe);
    R = sd^2*eye(length(Y));
    G = diag(1./lap_e2);

    C = zeros(size(meas_ind,1),size(fun2,3));
    for i=1:size(meas_ind,1)
        for j=1:size(fun2,3)
            C(i,j) = fun2(meas_ind(i,1),meas_ind(i,2),j);
        end
    end
    H = C*G;
    m_c  = P*H'/(H*P*H'+R)*Y
    mf_c = G*m_c;

    m_p  = zeros(size(uu));
    mf_p = zeros(size(uu));
    for i=1:length(m_c)
        m_p  = m_p + m_c(i) * fun2(:,:,i);
        mf_p = mf_p + mf_c(i) * fun2(:,:,i);
    end

    clf;
    subplot(2,2,1);
    pcolor(xx1,xx2,uu)
    shading flat;
    colormap('default')
    colorbar;
    title('Original input');
    cx = caxis;

    subplot(2,2,2);
    pcolor(xx1,xx2,m_p)
    shading flat;
    colormap('default')
    colorbar;
    title('Reconstructed input');
    caxis(cx);

    subplot(2,2,3);
    pcolor(xx1,xx2,ff_p)
    shading flat;
    colormap('default')
    colorbar;
    title('Original f(x)');
    cx = caxis;

    subplot(2,2,4);
    pcolor(xx1,xx2,mf_p)
    shading flat;
    colormap('default')
    colorbar;
    title('Reconstructed f(x)');
    caxis(cx);

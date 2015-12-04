%
% Simulate data from the PDE model (run demoPDEwithRBF first).
%

    F_sim = zeros(size(xgrid));
    U_sim = zeros(size(xgrid));

    m = zeros(size(F,1),1);
    P = blkdiag(0*eye(size(Fpd,1)),Pgp);
    C = blkdiag(0*eye(size(Fpd,1)),chol(Pgp,'lower'));
    
    x_sim = m + C * randn(size(F,1),1);
    C = chol(Q,'lower');
    
    XX = zeros(size(x_sim,1),length(t));
    
    for k=1:length(t)
        if k > 1
            x_sim = A * x_sim + C * randn(size(F,1),1);
        end
        
        XX(:,k) = x_sim;
        F_sim(:,k) = Hp * x_sim;
        U_sim(:,k) = Hup * x_sim;
    end

    sigma2 = 0.01*var(F_sim(:))
    yk_sim = zeros(size(yk));

    for k=1:length(t)
        % Evaluate the basis at the measurement points
        ind = find(abs(txk(:,1) - t(k)) < eps);
        meas_x = txk(ind,2);
    
        if length(ind) > 0
            E = zeros(length(meas_x),N);
            for i=1:length(meas_x)
                for j=1:N
                    xx = meas_x(i);
                    E(i,j) = eigenfun(j,xx);
                end
            end
            Hpd  = E;       % Projects into measurement points
            H = [Hpd zeros(size(Hpd,1),size(Fgp,1))];
            
            y = H * XX(:,k) + sqrt(sigma2) * randn(length(ind),1);
            yk_sim(ind) = y;
        end
    end

    f = F_sim(:);
    u = U_sim(:);
    yk = yk_sim;
    
    % First we plot the output
    clf;
    subplot(1,2,1);
    surf(tgrid, xgrid, F_sim)
    hold on
    plot3(txk(:,1), txk(:, 2), yk,  '.k', 'markersize', 15, 'linewidth', 1.5 )
    xlabel('Time')
    ylabel('Space')
    zlabel('Output function')
    title('Sample from the output function')

    % Now we plot the input
    subplot(1,2,2);
    surf(tgrid, xgrid, U_sim)    
    xlabel('Time')
    ylabel('Space')
    zlabel('input function')
    title('Sample for the input function')
    
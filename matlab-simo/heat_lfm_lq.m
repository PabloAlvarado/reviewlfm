%
% Controlled heat equation with basic LQ control. You should
% have run heat_basis_ex first and use the full basis set.
%

    %%
    % Random seed and parameters
    %
    s = RandStream('mt19937ar', 'Seed', 1e8);
    RandStream.setGlobalStream(s);

    %%
    % Design the LQ controller. This needs the estimation
    % to have the same size of basis as the simulation!
    %
    if size(Fpd,1) ~= size(F_sim)
        error('LQ cannot have different number of basis functions in simulation as estimation!');
    end
    this_should_be_zero = norm(Fpd-F_sim)
    
    X_lq = eye(size(Lpd,1));
    X_lq = blkdiag(X_lq,zeros(size(F,1)-size(X_lq,1)));
    U_lq = 0.05*eye(size(Lpd,2));
    [G,S_lq,E_lq] = lqr(F,Lc,X_lq,U_lq);
    
    %%
    % Solve the controlled equation
    %
    kf_x_mu_lfm_lq = zeros(size(xx1,1),size(xx1,2),length(tt));
    kf_u_mu_lfm_lq = zeros(size(xx1,1),size(xx1,2),length(tt));
    
    m = zeros(size(F,1),1);
    %m(1) = x0;
    %m(2) = dx0;
    P = blkdiag(1*eye(size(Fpd,1)),Pgp);
    P0 = P;

    MM = zeros(size(m,1),length(tt));
    PP = zeros(size(P,1),size(P,2),length(tt));
    
    plot_list = 1:10:length(tt);

    R = sd^2*eye(size(YYT,1));
    
    FFT_c_lfm_lq = zeros(size(UUT_c));
    FFT_p_lfm_lq = zeros(size(UUT_p));

    YYT_lfm_lq = zeros(size(meas_ind,1),length(tt));
    
    UUT_c_lfm_lq = zeros(size(UUT_c));
    UUT_p_lfm_lq = zeros(size(UUT_p));
    
    ff_c = zeros(size(FFT_c_lfm_lq,1),1);
    for k=1:size(FFT_c_lfm_lq,2)
        u = -G*m;
        
        UUT_c_lfm_lq(:,k) = u;
        
        for i=1:size(UUT_c_lfm_lq,1)
            UUT_p_lfm_lq(:,:,k) = UUT_p_lfm_lq(:,:,k) + UUT_c_lfm_lq(i,k) * fun(:,:,i);
        end
        
        
        % Actual system
        ff_c = A_sim*ff_c + B_sim*UUT_c(:,k) + B_sim*u;
        FFT_c_lfm_lq(:,k) = ff_c;
        for i=1:size(FFT_c_lfm_lq,1)
            FFT_p_lfm_lq(:,:,k) = FFT_p_lfm_lq(:,:,k) + FFT_c_lfm_lq(i,k) * fun(:,:,i);
        end
        for i=1:size(meas_ind,1)
            YYT_lfm_lq(i,k) = FFT_p_lfm_lq(meas_ind(i,1),meas_ind(i,2),k) + sd * randn;
        end
        
        % Estimate of the system
        m = A*m + Lc*u*dt;
        P = A*P*A' + Q;

        y = YYT_lfm_lq(:,k);
        S = H*P*H' + R;
        K = P * H' / S;
        m = m + K * (y - H*m);
        P = P - K * S * K';
        
        MM(:,k) = m;
        PP(:,:,k) = P;
        
        m_c = Hc*m;
        mu_c = Huc*m;
        
        for i=1:size(m_c,1)
            kf_x_mu_lfm_lq(:,:,k) = kf_x_mu_lfm_lq(:,:,k) + m_c(i) * fun2(:,:,i);
        end
        
        for i=1:size(mu_c,1)
            kf_u_mu_lfm_lq(:,:,k) = kf_u_mu_lfm_lq(:,:,k) + mu_c(i) * fun2(:,:,i);
        end
        
        if any(plot_list == k)
            subplot(2,2,1);
            hold off;
            pcolor(xx1,xx2,FFT_p_lfm_lq(:,:,k))
            caxis([0 mx]);
            shading flat;
            colormap('hot')
            colorbar;
            hold on;
            plot3(meas_x1,meas_x2,YYT_lfm_lq(:,k),'x');

            title(sprintf('Frame %d/%d',k,size(UUT,3)));

            subplot(2,2,2);
            pcolor(xx1,xx2,kf_x_mu_lfm_lq(:,:,k))
            caxis([0 mx]);
            shading flat;
            colormap('hot')
            colorbar;
            
            title('Kalman field');
            
            subplot(2,2,3);
            pcolor(xx1,xx2,UUT_p(:,:,k))
%            caxis([0 mx]);
            shading flat;
            colormap('hot')
            colorbar;
            
            title('Force input');
            
            subplot(2,2,4);
            pcolor(xx1,xx2,UUT_p_lfm_lq(:,:,k))
%            caxis([0 mx]);
            shading flat;
            colormap('hot')
            colorbar;
            
            title('Control signal');

            drawnow;

            %pause(0.01);
        end
    end

    %%
    % Plot the maxima (or means) at each time
    %
    maxim_r_lfm_lq = zeros(1,size(kf_x_mu_lfm_lq,3));
    maxim_e_lfm_lq = zeros(1,size(kf_x_mu_lfm_lq,3));
    for k=1:1:size(kf_x_mu_lfm_lq,3)
%        maxim_r_lfm_lq(k) = max(max(FF_p_lq(:,:,k)));
%        maxim_e_lfm_lq(k) = max(max(kf_x_mu_lq(:,:,k)));
        maxim_r_lfm_lq(k) = mean(mean(FFT_p_lfm_lq(:,:,k)));
        maxim_e_lfm_lq(k) = mean(mean(kf_x_mu_lfm_lq(:,:,k)));
    end
    
    clf;
    plot(tt,maxim_r_lfm_lq,tt,maxim_e_lfm_lq)
    
    %%
    % Compare to basic LQ
    %
    clf;
    plot(tt,maxim_r_lq,'--',tt,maxim_e_lq,'--',...
         tt,maxim_r_lfm_lq,tt,maxim_e_lfm_lq);
    legend('f LQ sim','f LQ kf','f LFM-LQ sim','f LFM-LQ kf');
    grid on;
    
    %%
    % Show a typical control signal
    %
    k = round(length(tt)/3);
    pcolor(xx1,xx2,UUT_p_lfm_lq(:,:,k))
    shading flat;
    colormap('hot')
    colorbar;
    title(sprintf('Control signal at time %.1f',tt(k)));
    
    
    
    
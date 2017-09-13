%
% Generate figures and stuff for the heat example
%

    %%
    % Form the model
    %
    heat_basis_ex;

    kf_x_mu0 = kf_x_mu;
    kf_u_mu0 = kf_u_mu;
    
    %%
    % Basic LQ
    %
    heat_basic_lq;
    
    
    %%
    % LFM LQ
    %
    heat_lfm_lq;
    
    %%
    % Plots
    %
    k = 70;
    tt(k)
    
    clf
    pcolor(xx1,xx2,FFT_p(:,:,k))
    caxis([0 mx]);
    shading flat;
    colormap('hot')
    colorbar;
    hold on;
    plot3(meas_x1,meas_x2,YYT(:,k),'o');

    print -depsc heat_open_x;
    
    %%
    clf
    pcolor(xx1,xx2,UUT_p(:,:,k))
    caxis([0 mx]);
    shading flat;
    colormap('hot')
    colorbar;
    
    print -depsc heat_open_u;

    %%
    clf
    pcolor(xx1,xx2,FFT_p_lq(:,:,k))
    caxis([0 mx]);
    shading flat;
    colormap('hot')
    colorbar;
    
    print -depsc heat_lq_x;
    
    %%
    clf
    pcolor(xx1,xx2,FFT_p_lfm_lq(:,:,k))
    caxis([0 mx]);
    shading flat;
    colormap('hot')
    colorbar;
%    hold on;
%    plot3(meas_x1,meas_x2,YYT(:,k),'o');

    print -depsc heat_lfm_lq_x;

    %%
    clf
    pcolor(xx1,xx2,UUT_p_lfm_lq(:,:,k))
    caxis([-mx 0]);
    shading flat;
    colormap('jet')
    colorbar;

    print -depsc heat_lfm_lq_c;
    
    %%
    clf;
    h = plot(tt,maxim_r_lq,'--',...
             tt,maxim_r_lfm_lq);
    legend('Basic LQR','LFM LQR');
    grid on;
    set(h(1),'LineWidth',2);
    set(h(1),'Color',[0.7 0.7 0.7]);
    set(h(2),'LineWidth',1);
    set(h(2),'Color',[0.0 0.0 0.0]);

    xlabel('Time{\it t}');
    ylabel('Maximum temperature{\it T(t)}');
 
    print -depsc heat_maxtemp;
    
    
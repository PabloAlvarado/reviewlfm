%
% Simulate data from heat equation in 2d using basis
% function approach.
%

    %%
    % Random seed and parameters
    %
    s = RandStream('mt19937ar', 'Seed', 1e8);
    RandStream.setGlobalStream(s);

    decay = 0.2;  % Parameters \lambda in the equation
    diffusion = 1e-3; % Diffusion constant
    sensitivity = 1; % Sensitivity parameter
    pos1 = [0.5 0.5]; % Heat source position start
    pos2 = [-1.5 -1.5]; % Heat source position end
    amp = 1;         % Heat source amplitude
    sd = 0.01;  % Measurement noise sd

    %%
    % Form the basis and evaluate it on a grid
    %

    m = 100; % Number of basis functions

    fprintf('Forming the basis of size %d ... ',m);
    LL = [1 1];
    [eigenfun,eigenval,NN] = domain_cartesian(m,2,LL);

    dd = [0.01 0.01];
    [xx1,xx2] = meshgrid((-LL(1)):dd(1):LL(1),(-LL(2)):dd(2):LL(2));

    fun = zeros(size(xx1,1),size(xx1,2),size(NN,1));

    for k=1:size(NN,1)
        nn = NN(k,:);
        fun(:,:,k) = reshape(eigenfun(nn,[xx1(:) xx2(:)]),size(xx1));
%        for i=1:size(xx1,1)
%            for j=1:size(xx1,2)
%                fun(i,j,k) = reshape(eigenfun(nn,[xx1(:) xx2(:)]),size(xx1));
%            end
%        end
    end
    
    lap_e = eigenval(NN);
    
    clf;
    pcolor(xx1,xx2,fun(:,:,end))
    shading flat;
    colormap('default')
    colorbar;

    pause(0.5)
    fprintf('Done.\n');

    %%
    % Check orthogonality
    %
    myinner = @(f1,f2) dd(1)*dd(2)*sum(sum(f1.*f2));

    C = zeros(size(fun,3));
    for i=1:size(fun,3)
        for j=1:size(fun,3)
            C(i,j) = myinner(fun(:,:,i),fun(:,:,j));
        end
    end

    clf;
    pcolor(C)
    shading flat;
    colormap('default')
    colorbar;
    
    pause(0.5)

    %%
    % Form the input function and project it onto the basis
    %
    dt = 0.1;
    tt = 0:dt:15;
    
    a1_t = 1 ./ (1 + exp(5*(1-tt)));
    a2_t = 1 ./ (1 + exp(5*(tt-tt(end)+7)));
    a_t = amp * a1_t .* a2_t;

    clf;
    subplot(2,2,1);
    plot(tt,a_t)

    z = 0.2;
    
    UUT_c = zeros(size(fun,3),length(tt));
    UUT   = zeros(size(xx1,1),size(xx1,2),length(tt));
    
    for k=1:size(UUT_c,2)
        s = k/size(UUT_c,2);
        pos = (1-s)*pos1 + s*pos2;
        u_x = 1 ./ ((xx1-pos(1)).^2 + (xx2-pos(2)).^2 + z^2);
        UUT(:,:,k) = a_t(k) * u_x;
    end
    
    subplot(2,2,2);
    pcolor(xx1,xx2,u_x);
    colorbar;
    shading interp;
    
    UUT_p = zeros(size(UUT));
    for k=1:size(UUT_c,2)
        for i=1:size(UUT_c,1)
            UUT_c(i,k) = myinner(UUT(:,:,k),fun(:,:,i));
            UUT_p(:,:,k) = UUT_p(:,:,k) + UUT_c(i,k) * fun(:,:,i);
        end
    end

    ind = round(length(tt)/2);
    
    subplot(2,2,3);
    pcolor(xx1,xx2,UUT(:,:,ind))
    ca = caxis;
    shading flat;
    colormap('hot')
    colorbar;
    title('Original input');

    subplot(2,2,4);
    pcolor(xx1,xx2,UUT_p(:,:,ind))
    caxis(ca);
    shading flat;
    colormap('hot')
    colorbar;
    title('Projected input');
    
    pause(0.5)


    %%
    % Animate the whole input
    %
    clf;
    mx = max(UUT(:));
    for k=1:size(UUT,3)
        pcolor(xx1,xx2,UUT(:,:,k))
        caxis([0 mx]);
        shading flat;
        colormap('hot')
        colorbar;
        title(sprintf('Frame %d/%d',k,size(UUT,3)));
        drawnow;
        %pause(0.01);
    end
    
    %%
    % Form the state-space model, solve the equation,
    % and form the measurements
    %
    F_sim = diag(-decay - diffusion * eigenval(NN));
    L_sim = eye(size(F_sim));
    Z = zeros(size(F_sim));
    tmp = expm([F_sim L_sim; Z Z]*dt);
    n = size(F_sim,1);
    A_sim = tmp(1:n,1:n);
    B_sim = tmp(1:n,n+1:end);
    
    dn1  = round(size(xx1,1)/20);
    dn2  = round(size(xx1,2)/20);
    ind1 = floor(dn1/2):dn1:ceil(size(xx1,1)-dn1/2);
    ind2 = floor(dn2/2):dn2:ceil(size(xx1,2)-dn2/2);
    
    [m1,m2] = meshgrid(ind1,ind2);
    meas_ind = [m1(:) m2(:)];
    meas_x1 = xx1(ind1,ind2);
    meas_x1 = meas_x1(:);
    meas_x2 = xx2(ind1,ind2);
    meas_x2 = meas_x2(:);
    
    FFT_c = zeros(size(UUT_c));
    FFT_p = zeros(size(UUT_p));

    YYT = zeros(size(meas_ind,1),length(tt));
    
    ff_c = zeros(size(FFT_c,1),1);
    for k=1:size(FFT_c,2)
        ff_c = A_sim*ff_c + B_sim*UUT_c(:,k);
        FFT_c(:,k) = ff_c;
        for i=1:size(FFT_c,1)
            FFT_p(:,:,k) = FFT_p(:,:,k) + FFT_c(i,k) * fun(:,:,i);
        end
        for i=1:size(meas_ind,1)
            YYT(i,k) = FFT_p(meas_ind(i,1),meas_ind(i,2),k) + sd * randn;
        end
        
        hold off;
        pcolor(xx1,xx2,FFT_p(:,:,k))
        caxis([0 mx]);
        shading flat;
        colormap('hot')
        colorbar;
        hold on;
        plot3(meas_x1,meas_x2,YYT(:,k),'x');
        
        title(sprintf('Frame %d/%d',k,size(UUT,3)));
        drawnow;
        %pause(0.01);
    end
    
    %%
    % Plot in 3d
    %
    k = round(length(tt)/2);
    
    hold off;
    surfl(xx1,xx2,FFT_p(:,:,k))
    shading flat;
    colormap('hot')
    hold on;
    alpha(0.8);
    plot3(meas_x1,meas_x2,YYT(:,k),'x');
    camproj perspective;
    

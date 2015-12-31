%
% Simulate data for 2D Poisson equation using basis
% function approach.
%

    %%
    % Random seed
    %
    s = RandStream('mt19937ar', 'Seed', 1e8);
    RandStream.setGlobalStream(s);

    %%
    % Form the basis and evaluate it on a grid
    %

    m = 200; % Number of basis functions

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

    clf;
    pcolor(xx1,xx2,fun(:,:,10))
    shading flat;
    colormap('default')
    colorbar;

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


    %%
    % Form the input function and project it onto the basis
    %
    c1 = [-0.3 0.3];
    s1 = 0.02;
    c2 = [0.2 -0.4];
    s2 = 0.02;
    uu = 100*exp(-((xx1-c1(1)).^2+(xx2-c1(2)).^2)/2/s1) ...
       + 100*exp(-((xx1-c2(1)).^2+(xx2-c2(2)).^2)/2/s2);

    uu_c = zeros(size(fun,3),1);
    uu_p = zeros(size(uu));
    for i=1:length(uu_c)
        uu_c(i) = myinner(uu,fun(:,:,i));
        uu_p = uu_p + uu_c(i) * fun(:,:,i);
    end

    clf;
    subplot(1,2,1);
    pcolor(xx1,xx2,uu)
    shading flat;
    colormap('default')
    colorbar;
    title('Original input');

    subplot(1,2,2);
    pcolor(xx1,xx2,uu_p)
    shading flat;
    colormap('default')
    colorbar;
    title('Projected input');

    %%
    % Solve the Poisson equation -D f = u and
    % generate the measurements
    %
    sd = 0.1;
    ff_c = zeros(size(fun,3),1);
    lap_e = eigenval(NN);
    ff_c = diag(1 ./ lap_e) * uu_c;
    ff_p = zeros(size(uu));
    for i=1:length(ff_c)
        ff_p = ff_p + ff_c(i) * fun(:,:,i);
    end

    dn   = round(size(xx1,1)/20);
    ind1 = round(dn/2):dn:round(size(xx1,1)-dn/2);
    ind2 = round(dn/2):dn:round(size(xx2,1)-dn/2);
    meas_ind = [];
    meas_x = [];
    Y = [];
    for i=1:length(ind1)
        for j=1:length(ind2)
            meas_ind = [meas_ind; ind1(i) ind2(j)];
            meas_x   = [meas_x; xx1(ind1(i),ind2(j)) xx2(ind1(i),ind2(j))];
            Y = [Y; ff_p(ind1(i),ind2(j))+sd*randn];
        end
    end


    clf;
    surf(xx1,xx2,ff_p)
    shading flat;
    colormap('default')
    alpha(0.9);
    colorbar;
    hold on;
    plot3(meas_x(:,1),meas_x(:,2),Y,'rx');
    camproj perspective;
    title('Solution of Poisson equation');


%
% Batch state-space solution to a PDE. You should have run
% demoPDEwithRBF_ss.m just before calling this.

    %%
    % Form the necessary covariances (takes a while)
    %
    Kfy = zeros(length(f),length(yk));
    Kuy = zeros(length(u),length(yk));
    Kyy = zeros(length(yk),length(yk));
    
    HH = cell(1,length(t));

    tic;
    fprintf('Evaluating Hs...\n');
    for k=1:length(t)
        % Evaluate the basis at measurement locations
        % for each time
        ind = find(abs(txk(:,1) - t(k)) < eps);
        meas_x = txk(ind,2);
        y = yk(ind);
        
%        E = zeros(length(meas_x),N);
%        for i=1:length(meas_x)
%            xx = meas_x(i);
%            %for j=1:N
%            %    E(i,j) = eigenfun(j,xx);
%            %end
%        end
        E = eigenfun(1:N,meas_x);
        Hpd  = E;       % Projects into measurement points
        H = [Hpd zeros(size(Hpd,1),size(Fgp,1))];
        HH{k} = H;
    end
    toc

    fprintf('Evaluating the cross-covariances...\n');
    
    tic;
    for i=1:length(t)
        fprintf('Processing %d/%d...\n',i,length(t));
        
        H = HH{i};
        for j=1:length(t)
            H2 = HH{j};
        
            tmp_fy = Hp*ss_crosscov(t(i),t(j),F,L,Qc,P0)*H2';
            tmp_uy = Hup*ss_crosscov(t(i),t(j),F,L,Qc,P0)*H2';
            tmp_yy = H*ss_crosscov(t(i),t(j),F,L,Qc,P0)*H2';
            tmp_ff = Hp*ss_crosscov(t(i),t(j),F,L,Qc,P0)*Hp';
            tmp_uu = Hup*ss_crosscov(t(i),t(j),F,L,Qc,P0)*Hup';
            
            ind_f1 = find(abs(tgrid - t(i)) < eps);
            ind_f2 = find(abs(tgrid - t(j)) < eps);
            ind_y1 = find(abs(txk(:,1) - t(i)) < eps);
            ind_y2 = find(abs(txk(:,1) - t(j)) < eps);

            Kfy(ind_f1,ind_y2) = tmp_fy;
            Kuy(ind_f1,ind_y2) = tmp_uy;
            Kyy(ind_y1,ind_y2) = tmp_yy;
            Kff(ind_f1,ind_f2) = tmp_ff;
            Kuu(ind_f1,ind_f2) = tmp_uu;
        end
    end
    Kyy = Kyy + sigma2 * eye(size(Kyy,1));
    toc
    
    %%
    % Compute the GP solutions
    %
    fprintf('Computing GP solution.\n');
    tic
    gp_u_mu = Kuy / Kyy * yk;
    %gp_u_V  = diag(Kuu - Kuy / Kyy * Kuy');
    gp_x_mu = Kfy / Kyy * yk;
    %gp_x_V  = diag(Kff - Kfy / Kyy * Kfy');
    toc
     
    %%
    % Final plot
    %
    subplot(2,2,1);
    pcolor(xgrid,tgrid,reshape(gp_u_mu,size(xgrid)));
    shading interp;
    title('sGP u');
    colorbar;
    
    subplot(2,2,2);
    pcolor(xgrid,tgrid,UU - reshape(gp_u_mu,size(xgrid)));
    shading interp;
    title('sGP u err');
    colorbar;
    
    subplot(2,2,3);
    pcolor(xgrid,tgrid,reshape(mean_pred_u,size(xgrid)));
    shading interp;
    title('GP u');
    colorbar;
    
    subplot(2,2,4);
    pcolor(xgrid,tgrid,UU - reshape(mean_pred_u,size(xgrid)));
    shading interp;
    title('GP u err');
    colorbar;

    mse_error_gp = mean((gp_u_mu(:) - u).^2)
    
    ks_sgp_diff = max(abs(gp_u_mu(:) - ks_u_mu(:)))

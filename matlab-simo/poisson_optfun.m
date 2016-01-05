function nlog_lh = poisson_optfun(s_ell,fun2,Y,sd,lap_e2,meas_ind)

    s = exp(s_ell(1));
    ell = exp(s_ell(2));

%    se_cov = @(norm_x) s^2 * exp(-norm_x.^2/2/ell^2);
    se_spec = @(norm_w) s^2 * (2*pi) * ell^2 * exp(-ell^2 * norm_w.^2/2);

    Spe = se_spec(sqrt(lap_e2));

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
    S = (H*P*H'+R);
    m_c  = P*H'/S*Y;
    mf_c = G*m_c;

%    nlog_lh = 0.5*log(det(S)) + 0.5*Y'/S*Y;    
    nlog_lh = 0.5*2*sum(log(diag(chol(S)))) + 0.5*Y'/S*Y;    
end


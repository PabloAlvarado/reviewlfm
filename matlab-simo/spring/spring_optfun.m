function nlog_lh = spring_se_optfun(s_ell,sd,b,dt,Y)
% nlog_lh = spring_se_optfun(s_ell,sd,b,dt,Y)

    s = s_ell(1);
    ell = s_ell(2);

    [Bgp,Agp] = se_pade(4,8,s,ell);
    [Fgp,Lgp,qgp,Hgp] = ratspec_to_ss(Bgp,Agp);
    [Fgp,Lgp,Hgp] = ss_balance(Fgp,Lgp,Hgp);

    Pgp = lyapchol(Fgp,Lgp*sqrt(qgp));
    Pgp = Pgp' * Pgp;
    
    Fsp = [0 1; -1 -b];
    Lsp = [0;1];
    Hsp = [1 0];
    
    Fjm = blkdiag(Fsp,Fgp);
    d = size(Fsp,1);
    Fjm(1:d,d+1:end) = Lsp * Hgp;
    Ljm  = [zeros(d,1); Lgp];
    Hjm  = [Hsp zeros(1,size(Fgp,1))];
    qjm = qgp;
    R = sd^2;

    m = zeros(size(Fjm,1),1);
    P = blkdiag(eye(size(Fsp,1)),Pgp);
    
    [Ajm,Qjm] = lti_disc(Fjm,Ljm,qjm,dt);
    
    nlog_lh = 0;
    
    for k=1:length(Y)
        m = Ajm*m;
        P = Ajm*P*Ajm' + Qjm;
        
        if ~isnan(Y(k))
            S = Hjm*P*Hjm' + R;
            K = P * Hjm' / S;
            v = (Y(k) - Hjm*m);
            m = m + K * v;
            P = P - K * S * K';

            nlog_lh = nlog_lh + 0.5*log(S) + 0.5*v^2/S;
        end
        
    end
    

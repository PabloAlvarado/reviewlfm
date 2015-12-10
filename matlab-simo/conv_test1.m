%
% Test the convergence of Taylor/Pade approximations
%

% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

% Copyright (C) 2014 Simo Sarkka

%%
% Parameters and exact covariance functions and spectra
%
    s = 1;
    ell = 1;

    se_spec = @(w) s^2 * sqrt(2*pi) * ell * exp(-ell^2 * w.^2/2);
    se_cov  = @(t) s^2 * exp(-t.^2/2/ell^2);

    ctay = ((pi-1)/sqrt(2))^2;
    
%%
% Compute errors for plotting
%

    % Turn off the almost-singularity warnings that we get.
    warning off MATLAB:nearlySingularMatrix;

    err_list1 = [];
    err_list2 = [];
    err_list3 = [];
    err_list4 = [];
    err_list5 = [];
    err_list6 = [];

    rmse_list1 = [];
    rmse_list2 = [];
    rmse_list3 = [];
    rmse_list4 = [];
    rmse_list5 = [];
    rmse_list6 = [];

    ord_list1 = [];
    ord_list2 = [];
    ord_list3 = [];
    ord_list4 = [];
    ord_list5 = [];
    ord_list6 = [];
    
    cond_list1 = [];
    cond_list2 = [];
    cond_list3 = [];
    cond_list4 = [];
    cond_list5 = [];
    cond_list6 = [];

    tau = -3:0.01:3;
    
    nmax = 20;
    for n=1:nmax
        [B,A] = se_taylor(n,s,ell);    
        [F,L,q,H] = ratspec_to_ss(B,A);
        [F,L,H] = ss_balance(F,L,H);
%        err = abs(se_cov(0) - ss_cov(0,F,L,q,H));
        err = max(abs(se_cov(tau) - ss_cov(tau,F,L,q,H)));
        rmse = sqrt(mean((se_cov(tau) - ss_cov(tau,F,L,q,H)).^2));
        rmse_list1(n) = rmse;
        err_list1(n) = err;
        ord_list1(n) = size(F,1);
        cond_list1(n) = cond(F);
        
        [B,A] = se_power_taylor(n,1,s,ell);    
        [F,L,q,H] = ratspec_to_ss(B,A);
        [F,L,H] = ss_balance(F,L,H);
%        err = abs(se_cov(0) - ss_cov(0,F,L,q,H));
        err = max(abs(se_cov(tau) - ss_cov(tau,F,L,q,H)));
        rmse = sqrt(mean((se_cov(tau) - ss_cov(tau,F,L,q,H)).^2));
        rmse_list2(n) = rmse;
        err_list2(n) = err;
        ord_list2(n) = size(F,1);
        cond_list2(n) = cond(F);
        
        [B,A] = se_power_taylor(n,2,s,ell);    
        [F,L,q,H] = ratspec_to_ss(B,A);
        if size(F,1) <= nmax
            [F,L,H] = ss_balance(F,L,H);
%            err = abs(se_cov(0) - ss_cov(0,F,L,q,H));
            err = max(abs(se_cov(tau) - ss_cov(tau,F,L,q,H)));
            rmse = sqrt(mean((se_cov(tau) - ss_cov(tau,F,L,q,H)).^2));
            rmse_list3(n) = rmse;
            err_list3(n) = err;
            ord_list3(n) = size(F,1);
            cond_list3(n) = cond(F);
        end
            
        [B,A] = se_power_pade(n,2,4,s,ell);
        [F,L,q,H] = ratspec_to_ss(B,A);
        if size(F,1) <= nmax
            [F,L,H,T] = ss_balance(F,L,H);
%            err = abs(se_cov(0) - ss_cov(0,F,L,q,H));
            err = max(abs(se_cov(tau) - ss_cov(tau,F,L,q,H)));
            rmse = sqrt(mean((se_cov(tau) - ss_cov(tau,F,L,q,H)).^2));
            rmse_list4(n) = rmse;
            err_list4(n) = err;
            ord_list4(n) = size(F,1);
            cond_list4(n) = cond(F);
        end
        
        [B,A] = se_pade(2*n,4*n,s,ell);
        [F,L,q,H] = ratspec_to_ss(B,A);
        if size(F,1) <= nmax
            [F,L,H,T] = ss_balance(F,L,H);
%            err = abs(se_cov(0) - ss_cov(0,F,L,q,H));
            err = max(abs(se_cov(tau) - ss_cov(tau,F,L,q,H)));
            rmse = sqrt(mean((se_cov(tau) - ss_cov(tau,F,L,q,H)).^2));
            rmse_list5(n) = rmse;
            err_list5(n) = err;
            ord_list5(n) = size(F,1);
            cond_list5(n) = cond(F);
        end
        
        [B,A] = se_power_mtaylor(n,s,ell,ctay);
        [F,L,q,H] = ratspec_to_ss(B,A);
        if size(F,1) <= nmax
            [F,L,H] = ss_balance(F,L,H);
%            err = abs(se_cov(0) - ss_cov(0,F,L,q,H));
            err = max(abs(se_cov(tau) - ss_cov(tau,F,L,q,H)));
            rmse = sqrt(mean((se_cov(tau) - ss_cov(tau,F,L,q,H)).^2));
            rmse_list6(n) = rmse;
            err_list6(n) = err;
            ord_list6(n) = size(F,1);
            cond_list6(n) = cond(F);
        end
    end
    
    fprintf('Done.\n');

%%  
% Plot the final figure
%
    clf;
    set(gcf,'PaperType','a4');
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperPosition',[0.25 2.5 4 3]);
    set(gca,'Position',[0.25 2.5 4 3]);

    clf;
    h = semilogy(ord_list1,err_list1,'*-',...
             ord_list2,err_list2,'+-',...
             ord_list3,err_list3,'^-',...
             ord_list4,err_list4,'v-',...
             ord_list5,err_list5,'x-',...
             ord_list6,err_list6,'o-');
    ylim([10^-9 1]);
    set(h,'Markersize',10);
    set(h,'Linewidth',2);
    legend('Taylor n','Taylor [1]^n','Taylor [2]^n','Pade [2/4]^n','Pade [2n/4n]','MTaylor^n',3);
    grid on;
    xlabel('State dimension');
    ylabel('Maximum covariance function error');
    
    print -depsc cov_max;


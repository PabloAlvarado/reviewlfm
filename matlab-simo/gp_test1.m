%
% Test state-space GP regression with SE covariance
%

% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

% Copyright (C) 2014 Simo Sarkka

%%
% Create the data
%

randn('state',0);
rand('state',0);

dt=1;
dtt=0.01;
sd=0.1;
TT = (1:1100) * dtt;
mi = 100:100:1000;
T = TT(mi);
X  = sin(T);
XX = sin(TT);
Y  = X + sd * randn(size(X));

clf;
plot(TT,XX,'--',T,Y,'o');
grid on;
axis([0 10 -1.5 1.5]);

%%
% Parameters
%
s = 1;
ell = 1;

se_spec = @(w) s^2 * sqrt(2*pi) * ell * exp(-ell^2 * w.^2/2);
se_cov  = @(t) s^2 * exp(-t.^2/2/ell^2);


%%
% Do basic GP regression
%

        k = length(Y);
        K11 = zeros(k,k);
        for i=1:k
            for j=1:k
                K11(i,j) = se_cov(T(i)-T(j));
            end
        end
        K21 = zeros(length(TT),k);
        for i=1:size(K21,1)
            for j=1:k
                K21(i,j) = se_cov(TT(i) - T(j));
            end
        end
        mu = (K21 / (K11 + eye(k)*sd^2) * Y(1:k)')';
        V  = se_cov(0) - diag(K21 / (K11 + eye(k)*sd^2) * K21')';

%%
% Plot the above
%
    clf;
    set(gcf,'PaperType','a4');
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperPosition',[0.25 2.5 4 3]);
    set(gca,'Position',[0.25 2.5 4 3]);

    clf;
    p=patch([TT fliplr(TT)], ...
            [mu + 1.96*sqrt(V) fliplr(mu - 1.96*sqrt(V))],1);
    set(p,'EdgeColor','none')
    colormap(0.8*[1 1 1]);
    hold on;
    h = plot(TT,XX,'k--',TT,mu,'b-',T,Y,'ro');
    set(h,'Linewidth',3);
    set(h(1),'Linewidth',2);
    set(h,'MarkerSize',10);
    grid on;
    axis([0 11 -2.5 2]);
    legend(h,'True function','Regression mean','Observations',3);

    print -depsc reg_ex;


%%
% Test SEs in state-space form
%

% Turn off the almost-singularity warnings that we get.
warning off MATLAB:nearlySingularMatrix;

order_list = 1:20;

ctay = ((pi-1)/sqrt(2))^2;

se_fun_all = {...
    @(order) se_taylor(order,s,ell),...
    @(order) se_power_taylor(order,1,s,ell),...
    @(order) se_power_taylor(order,2,s,ell),...
    @(order) se_power_pade(order,2,4,s,ell),...
    @(order) se_pade(2*order,4*order,s,ell),...
    @(order) se_power_mtaylor(order,s,ell,ctay),...
};

ks_mu_all = cell(length(order_list),length(se_fun_all));
ks_V_all  = cell(length(order_list),length(se_fun_all));;
dim_list_all = {};
mean_err_all = {};
var_err_all  = {};

for n=1:length(se_fun_all)
    fprintf('Running %d/%d\n',n,length(se_fun_all));
    
    se_fun = se_fun_all{n};
    dim_list = [];
    mean_err = [];
    var_err  = [];

    for order=order_list
        
        [B,A] = se_fun(order);
        [F,L,q,H] = ratspec_to_ss(B,A);
        if size(F,1) <= order_list(end)
            [F,L,H] = ss_balance(F,L,H);
            Pinf = lyapchol(F,L*sqrt(q));
            Pinf = Pinf' * Pinf;
        
            kf_mu = [];
            kf_V  = [];
        
            m = zeros(size(F,1),1);
            P = Pinf;
        
            MM = zeros(size(m,1),length(TT));
            PP = zeros(size(P,1),size(P,2),length(TT));
        
            [A,Q] = lti_disc(F,L,q,TT(2)-TT(1));
        
            count = 0;
            for k=1:length(TT)
                m = A*m;
                P = A*P*A' + Q;
                i = find(mi == k);
                
                if ~isempty(i)
                    y = Y(i);
                    S = H*P*H' + sd^2;
                    K = P * H' / S;
                    m = m + K * (y - H*m);
                    P = P - K * S * K';
                    
                end
                
                MM(:,k) = m;
                PP(:,:,k) = P;
                
                kf_mu = [kf_mu H*m];
                kf_V  = [kf_V  H*P*H'];
                
            end
            
            MS = MM;
            PS = PP;
            ms = m;
            Ps = P;
            
            ks_mu = H*m;
            ks_V  = H*P*H';
            
            for k=size(MM,2)-1:-1:1
                m = MM(:,k);
                P = PP(:,:,k);
                
                mp = A*m;
                Pp = A*P*A' + Q;
                
                G = P*A'/Pp;
                ms = m + G * (ms - mp);
                Ps = P + G * (Ps - Pp) * G';
                
                MS(:,k) = ms;
                PS(:,:,k) = Ps;
                
                ks_mu = [H*ms ks_mu];
                ks_V  = [H*Ps*H' ks_V];
                
            end
            
            err = max(abs(ks_mu - mu));
            mean_err = [mean_err err];
            
            err = max(abs(ks_V - V));
            var_err  = [var_err err];
            
            dim_list = [dim_list length(m)];
            
            ks_mu_all{n,order} = ks_mu;
            ks_V_all{n,order}  = ks_V;
        end 
    end
    
    dim_list_all{n} = dim_list;
    mean_err_all{n} = mean_err;
    var_err_all{n}  = var_err;
end
fprintf('Done.\n');

%%
% Final plot for mean
%

    clf;
    set(gcf,'PaperType','a4');
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperPosition',[0.25 2.5 4 3]);
    set(gca,'Position',[0.25 2.5 4 3]);

    clf;
    h = semilogy(dim_list_all{1},mean_err_all{1},'*-',...
             dim_list_all{2},mean_err_all{2},'+-',...
             dim_list_all{3},mean_err_all{3},'^-',...
             dim_list_all{4},mean_err_all{4},'v-',...
             dim_list_all{5},mean_err_all{5},'x-',...
             dim_list_all{6},mean_err_all{6},'o-');
    set(h,'Markersize',10);
    set(h,'Linewidth',2);
    ylim([10^-9 1]);
    legend('Taylor n','Taylor [1]^n','Taylor [2]^n','Pade [2/4]^n','Pade [2n/4n]','MTaylor^n',3);
    grid on;
    xlabel('State dimension');
    ylabel('Maximum mean error');

    print -depsc mean_max;

%%
% Final plot for variance
%
    clf;
    set(gcf,'PaperType','a4');
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperPosition',[0.25 2.5 4 3]);
    set(gca,'Position',[0.25 2.5 4 3]);
    
    clf;
    h = semilogy(dim_list_all{1},var_err_all{1},'*-',...
             dim_list_all{2},var_err_all{2},'+-',...
             dim_list_all{3},var_err_all{3},'^-',...
             dim_list_all{4},var_err_all{4},'v-',...
             dim_list_all{5},var_err_all{5},'x-',...
             dim_list_all{6},var_err_all{6},'o-');
    set(h,'Markersize',10);
    set(h,'Linewidth',2);
    ylim([10^-9 1]);
    legend('Taylor n','Taylor [1]^n','Taylor [2]^n','Pade [2/4]^n','Pade [2n/4n]','MTaylor^n',3);
    grid on;
    xlabel('State dimension');
    ylabel('Maximum variance error');

    print -depsc var_max;


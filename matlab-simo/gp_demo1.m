%
% Demonstrate GP regression with Matern covariance
%

%%
% Create the data
%

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
% Do basic GP regression
%
magnSigma2 = 1;
lengthScale = 1;
matern = @(tau) magnSigma2 * ...
    exp(-sqrt(3)*abs(tau)./lengthScale).*(sqrt(3)*abs(tau)/lengthScale+1);

count = 0;

clf;
for k=0:length(Y)
    if k == 0
        mu = zeros(size(TT));
        V = matern(0) * ones(size(TT));
    else
        K11 = zeros(k,k);
        for i=1:k
            for j=1:k
                K11(i,j) = matern(T(i)-T(j));
            end
        end
        K21 = zeros(length(TT),k);
        for i=1:size(K21,1)
            for j=1:k
                K21(i,j) = matern(TT(i) - T(j));
            end
        end
        mu = (K21 / (K11 + eye(k)*sd^2) * Y(1:k)')';
        V  = matern(0) - diag(K21 / (K11 + eye(k)*sd^2) * K21')';
    end

    clf;
    p=patch([TT fliplr(TT)], ...
            [mu + 1.96*sqrt(V) fliplr(mu - 1.96*sqrt(V))],1);
    set(p,'EdgeColor','none')
    hold on;
    alpha(0.5);
    h = plot(TT,XX,'k--',T,Y,'ro',TT,mu,'k-');
    set(h,'Linewidth',2);
    set(h,'MarkerSize',10);
    grid on;
    axis([0 11 -1.5 1.5]);

end

plot(TT,mu,TT,mu-1.96*sqrt(V),'--',TT,mu+1.96*sqrt(V),'--');

%%
% Matern in state-space form
%
    % Derived constants
    lambda = sqrt(3)/lengthScale;
  
    % Feedback matrix
    F = [0,          1;
         -lambda^2,  -2*lambda];

    % Noise effect matrix
    L = [0;   1];

    % Spectral density
    Qc = 12*sqrt(3)/lengthScale^3*magnSigma2;

    % Observation model
    H = [1,   0];
  
    Pinf = are(F',zeros(size(F)),L*Qc*L');

    kf_mu = [];
    kf_V  = [];
    
    m = [0;0];
    P = Pinf;

    MM = zeros(size(m,1),length(TT));
    PP = zeros(size(P,1),size(P,2),length(TT));
    
    [A,Q] = lti_disc(F,L,Qc,TT(2)-TT(1));
    
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
        
        kf_mu = [kf_mu m(1)];
        kf_V  = [kf_V  P(1,1)];
        
        if ~isempty(i) || k == length(TT)
            clf;
            p=patch([TT(1:k) fliplr(TT(1:k))], ...
                    [kf_mu + 1.96*sqrt(kf_V) fliplr(kf_mu - 1.96*sqrt(kf_V))],1);
            set(p,'EdgeColor','none')
            hold on;
            alpha(0.5);
            h = plot(TT,XX,'k--',T,Y,'ro',TT(1:k),kf_mu,'k-');
            set(h,'Linewidth',2);
            set(h,'MarkerSize',10);
            grid on;
            axis([0 11 -1.5 1.5]);
            
        end
    end
    
    MS = MM;
    PS = PP;
    ms = m;
    Ps = P;
    
    ks_mu = m(1);
    ks_V  = P(1,1);

    count = 0;
    
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
        
        ks_mu = [ms(1) ks_mu];
        ks_V  = [Ps(1,1) ks_V];
        
        i = find(mi == k);
        if ~isempty(i) || k == 1
            clf;
            p=patch([TT(k:end) fliplr(TT(k:end))], ...
                    [ks_mu + 1.96*sqrt(ks_V) fliplr(ks_mu - 1.96*sqrt(ks_V))],1);
            set(p,'EdgeColor','none')
            hold on;
            alpha(0.5);
            h = plot(TT,XX,'k--',T,Y,'ro',TT(k:end),ks_mu,'k-');
            set(h,'Linewidth',2);
            set(h,'MarkerSize',10);
            grid on;
            axis([0 11 -1.5 1.5]);
            
        end
        
    end
    
    plot(TT,mu,TT,ks_mu-1.96*sqrt(V),'--',TT,ks_mu+1.96*sqrt(V),'--');

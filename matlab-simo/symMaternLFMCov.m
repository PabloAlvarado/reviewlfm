%
% Compute symbolic covariance function for a spring
% model with Matern 3/2 input via matrix exponential
% based methods.
%

    %%
    % Parameters. We have
    %
    %  k(t) = s^2 exp(-sqrt(3) |t| / ell) [sqrt(3) |t| / ell + 1]
    %
    % which corresponds to
    %
    %   dX/dt = [0 1; -lambda^2 -2*lambda] X + [0 1] W
    %       y = [1 0]
    %
    syms s ell positive;
    syms tau t tp real;
    
    magnSigma2 = s^2;
    lengthScale = ell;

    matern = @(tau) magnSigma2 * ...
    exp(-sqrt(3)*abs(tau)./lengthScale).*(sqrt(3)*abs(tau)/lengthScale+1);
    
    % Derived constants
    lambda = sqrt(3)/lengthScale;
  
    % Feedback matrix
    Fgp = [0,          1;
           -lambda^2,  -2*lambda]

    % Noise effect matrix
    Lgp = [0; 1]

    % Spectral density
    qgp = 12*sqrt(3)/lengthScale^3*magnSigma2

    % Observation model
    Hgp = [1 0]
    
    % Stationary covariance (in closed form)
    Pgp = diag([limit(matern(tau)) -limit(diff(matern(tau),2))])
    Fgp*Pgp+Pgp*Fgp' + Lgp*qgp*Lgp'
    
    %%
    % State space representation of the spring model:
    %
    %   x'' + lam x' + gam x = sen u
    %
    %  [x;x']' = [0 1; -gam -lam] + [0;sen] u
    %  y = [1 0] [x;x']
    %
    syms dc sc se positive;
    lam = dc;
    gam = sc;
    sen = se;
    
    Fsp = [0 1; -gam -lam]
    Lsp = [0;sen]
    Hsp = [1 0]

    %%
    % Create the joint state-space
    %
    F = blkdiag(Fsp,Fgp);
    d = size(Fsp,1);
    F(1:d,d+1:end) = Lsp * Hgp
    L = [zeros(d,1); Lgp]
    H = [Hsp zeros(1,size(Fgp,1))]
    Hu = [zeros(1,d) Hgp]
    q = qgp;
    
    %%
    % Consider the state-space model (i.e., SDE)
    %
    %   dx = F x dt + L dB, E[dB^2] = qc dt
    %
    % with some x(0) ~ N(0,P0). This how has the solution
    %
    %   x(t) = exp(F t) x(0)
    %   + int_0^t exp(F (t-s)) L dB(s)
    %
    % and hence
    %
    %   E[x(t) x^T(t + tau)]
    %   = exp(F t) P0 exp(F (t + tau))^T
    %   + int_0^min{t,t+tau} exp(F (t-s)) L qc L^T exp(F (t-s+tau))^T ds
    %
    % If we assume that tau > 0 then we have
    %
    %   E[x(t) x^T(t + tau)]
    %   = exp(F t) P0 exp(F t) exp(F tau)^T
    %   + int_0^t exp(F (t-s)) L qc L^T exp(F (t-s)) ds exp(F tau)^T
    %   = [A(t) P0 A^T(t) + Q(t)] exp(F tau)^T
    %
    % where A and Q are returned by the lti_disc.
    % Similarly we get (still for tau > 0):
    %
    %   E[x(t + tau) x^T(t)]
    %   = exp(F tau) exp(F t) P0 exp(F t) 
    %   + exp(F tau) int_0^t exp(F (t-s)) L qc L^T exp(F (t-s)) ds
    %   = exp(F tau) [A(t) P0 A^T(t) + Q(t)] 
    %
    A = simplify(ilaplace(inv(s*eye(size(F,1))-F)))
    
    %%
    n   = size(F,1);
    Phi = [F L*q*L'; zeros(n,n) -F']
    Psi = simplify(ilaplace(inv(s*eye(size(Phi,1))-Phi)))
    AB  = Psi*[zeros(n,n);eye(n)];
    Q   = AB(1:n,:)/AB((n+1):(2*n),:)
    
    
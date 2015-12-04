% C = ss_crosscov(t,s,F,L,q,P0)
%
% Evaluate cross-covariance function C(t,s) of state-space model.
% If P0 is not given, the stationary covariance is returned.

%
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
% The stationary stage is obtained by P0 = Pinf, which
% then gives A(t) P0 A^T(t) + Q(t) = Pinf

function C = ss_crosscov(t,s,F,L,q,P0)
    if nargin < 6
        P0 = [];
    end
    
    if s < t
        C = ss_crosscov(s,t,F,L,q,P0)';
    else
        if isempty(P0)
            Pinf = lyapchol(F,L*sqrt(q));
            Pinf = Pinf' * Pinf;
            C = Pinf*expm(F*(s-t))';
        else
            %[At,Qt] = lti_disc(F,L,q,t);
            % This uses the result of Davison (1968), requires
            % that a stable solution exists.
            At = expm(F*t);
            P = lyap(F,L*q*L');
            Qt = P - At*P*At';
            C = (At*P0*At' + Qt)*expm(F*(s-t))';
        end
    end
    
        
 
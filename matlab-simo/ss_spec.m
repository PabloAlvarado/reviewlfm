% S_approx = ss_spec(w,F,L,q,H)
%
% Spectral density of state space model.
function S_approx = ss_spec(w,F,L,q,H)

    S_approx = arrayfun(@(wt) H/(F-1i*wt*eye(length(F)))*L*q*L'/((F+1i*wt*eye(length(F))).')*H',w);
    S_approx = real(S_approx);
end


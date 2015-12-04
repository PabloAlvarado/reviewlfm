% cov_approx = ss_cov(tau,F,L,q,H)
%
% Covariance function of state space model.
function cov_approx = ss_cov(tau,F,L,q,H)

    %Pinf = are(F',zeros(size(F)),L*q*L');
%    Pinf = lyap(F,F',L*q*L');
    Pinf = lyapchol(F,L*sqrt(q));
    Pinf = Pinf' * Pinf;

    % Initialize covariance
    cov_approx = zeros(size(tau));
  
    % Evaluate positive parts
    cov_approx(tau >= 0) = arrayfun(@(taut) H*Pinf*expm(taut*F)'*H',tau(tau >= 0));
  
    % Evaluate negative parts
    cov_approx(tau < 0) = arrayfun(@(taut) H*expm(-taut*F)*Pinf*H',tau(tau < 0));
end


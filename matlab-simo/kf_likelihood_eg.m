function [edata,gdata] = kf_likelihood_eg(F,L,Qc,H,R,Pinf,dF,dQc,dR,dPinf,y,x)
%% kf_likelihood_eg - Log marginal likelihood and gradient
%
% Syntax:
%   [e,eg] = kf_likelihood_eg(F,L,Qc,H,R,Pinf,dF,dQc,dR,dPinf,y,x)
%
% In:
%   F        - Feedback matrix
%   L        - Noise effect matrix
%   Qc       - Spectral density of the white noise process
%   H        - Observation model matrix
%   R        - Measurement noise covariance
%   Pinf     - Covariance of the stationary process
%   dF       - Partial derivatives of F w.r.t. all parameters
%   dQc      - Partial derivatives of Qc w.r.t. all parameters
%   dR       - Partial derivatives of R w.r.t. all parameters
%   dPinf    - Partial derivatives of Pinf w.r.t. all parameters
%   y        - The observed data
%   x        - Observation time points
%
% Out:
%   edata    - Return the evaluated energy (negative marginal log-likelihood)
%   gdata    - The energy gradient
%   
% Description:
%   This function performs sequential inference using the Kalman filter
%   algorithm and returns the log marginal likelihood of the model given
%   the parameters theta. The state space model is
%
%     df(t)/dt = F(theta)*f(t) + L*w(t),
%     y_k      = H*f(t_k)      + r,       r ~ N(0,R(theta)),
%
%   where w(t) is temporally white noise with spectral density Qc(theta).
%   Furthermore, the likelihood gradient is also returned.
%
% Copyright: 
%   2013-2014 - Arno Solin (initial version)
%
% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  % Sort values
  [x,ind] = sort(x(:));
  y = y(ind);
        
  % State dimension, number of data points and number of parameters
  n            = size(F,1); 
  steps        = numel(y);
  nparam       = size(dF,3);
  jitterSigma2 = 1e-9;

  % Evaluate likelihood only
  if nargout==1

    % Initialize likelihood
    edata = 0;
    
    % Set up
    m  = zeros(n,1);
    P  = Pinf;
    dt = -inf;    
    
    % Run filter for evaluating the marginal likelihood
    for k=1:steps
        
      % The previous time step
      dt_old = dt;
        
      % The time discretization step length
      if (k>1)
        dt = x(k)-x(k-1);
      else
        dt = 0;
      end
      
      % Solve A using the method by Davison
      if abs(dt-dt_old) > 1e-9
        
        % Discrete-time solution (only for stable systems)
        A  = expm(F*dt);
        Q  = Pinf - A*Pinf*A';
                        
      end
      
      % Prediction step
      m = A * m;
      P = A * P * A' + Q;

      % Update step
      S = H*P*H'+R;
      K = P*H'/S;
      v = y(k)-H*m;
      m = m + K*v;
      P = P - K*H*P;
      
      % Update log likelihood
      edata = edata + 1/2*log(det(2*pi*S)) + 1/2*v'/S*v;
      
    end
   
  % Also evaluate the likelihood gradient
  else
    
    % Allocate space for results
    edata = 0;
    gdata = zeros(1,nparam);
    
    % Set up
    Z  = zeros(n);
    m  = zeros(n,1);
    P  = Pinf;
    dm = zeros(n,nparam);
    dP = dPinf;
    dt = -inf;
    
    % Allocate space for expm results
    AA = zeros(2*n,2*n,nparam);
    
    % Loop over all observations
    for k=1:steps
        
        % The previous time step
        dt_old = dt;
        
        % The time discretization step length
        if (k>1)
            dt = x(k)-x(k-1);
        else
            dt = 0;
        end
        
        % Loop through all parameters (Kalman filter prediction step)
        for j=1:nparam
            
            % Should we recalculate the matrix exponential?
            if abs(dt-dt_old) > 1e-9
                
                % The first matrix for the matrix factor decomposition
                FF = [ F        Z;
                      dF(:,:,j) F];
                
                % Solve the matrix exponential
                AA(:,:,j) = expm(FF*dt);
                
            end
            
            % Solve the differential equation
            foo     = AA(:,:,j)*[m; dm(:,j)];
            mm      = foo(1:n,:);
            dm(:,j) = foo(n+(1:n),:);
            
            % The discrete-time dynamical model
            if (j==1)
                A  = AA(1:n,1:n,j);
                Q  = Pinf - A*Pinf*A';
                PP = A*P*A' + Q;
            end
            
            % The derivatives of A and Q
            dA = AA(n+1:end,1:n,j);
            %dQ = dPinf(:,:,j) - dA*Pinf*A' - A*dPinf(:,:,j)*A' - A*Pinf*dA';
            dAPinfAt = dA*Pinf*A';
            dQ = dPinf(:,:,j) - dAPinfAt - A*dPinf(:,:,j)*A' - dAPinfAt';
            
            % The derivatives of P
            %dP(:,:,j) = dA*P*A' + A*dP(:,:,j)*A' + A*P*dA' + dQ;
            dAPAt = dA*P*A';
            dP(:,:,j) = dAPAt + A*dP(:,:,j)*A' + dAPAt' + dQ;
        end
        
        % Set predicted m and P
        m = mm;
        P = PP;
        
        % Start the Kalman filter update step and precalculate variables
        S = H*P*H' + R;
        [LS,notposdef] = chol(S,'lower');
        
        % If matrix is not positive definite, add jitter
        if notposdef>0,
            jitter = jitterSigma2*diag(rand(size(S,1),1));
            [LS,notposdef] = chol(S+jitter,'lower');
            
            % Return nan to the optimizer
            if notposdef>0,
                edata = nan; gdata = nan*gdata;
                return;
            end
        end
        
        % Continue update
        HtiS = H'/LS/LS';
        iS   = eye(size(S))/LS/LS';
        K    = P*HtiS;
        v    = y(k) - H*m;
        vtiS = v'/LS/LS';
        
        % Loop through all parameters (Kalman filter update step derivative)
        for j=1:nparam
            
            % Innovation covariance derivative
            dS = H*dP(:,:,j)*H' + dR(:,:,j);
            
            % Evaluate the energy derivative for j (optimized from above)
            gdata(j) = gdata(j) ...
                + .5*sum(iS(:).*dS(:)) ...
                - .5*(H*dm(:,j))*vtiS' ...
                - .5*vtiS*dS*vtiS'     ...
                - .5*vtiS*(H*dm(:,j));
            
            % Kalman filter update step derivatives
            dK        = dP(:,:,j)*HtiS - P*HtiS*dS/LS/LS';
            dm(:,j)   = dm(:,j) + dK*v - K*H*dm(:,j);
            dKSKt     = dK*S*K';
            dP(:,:,j) = dP(:,:,j) - dKSKt - K*dS*K' - dKSKt';
            
        end
        
        % Evaluate the energy
        edata = edata + .5*size(S,1)*log(2*pi) + sum(log(diag(LS))) + .5*vtiS*v;
        %edata = edata + 1/2*log(det(2*pi*S)) + 1/2*v'/S*v;
        
        % Finish Kalman filter update step
        m = m + K*v;
        P = P - K*S*K';
        
    end
    
  end
    
  
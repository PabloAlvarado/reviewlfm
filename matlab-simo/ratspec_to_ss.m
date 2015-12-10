% [F,L,q,H,Pinf] = ratspec_to_ss(B,A,rtype)
%
% rtype = 0 means controllable canonical form and
% rtype = 1 the observable one.

% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

% Copyright (C) 2014 Simo Sarkka

function [F,L,q,H,Pinf] = ratspec_to_ss(B,A,rtype)
    
    if nargin < 3
        rtype = 0;
    end

    q = polyval(B,0) / polyval(A,0);
    
    LB = B ./ 1i.^(length(B)-1:-1:0);
    LA = A ./ 1i.^(length(A)-1:-1:0);

    BR = roots(LB);
    AR = roots(LA);
    
    GB = poly(BR(real(BR) < 0));
    GA = poly(AR(real(AR) < 0));
    
    GB = GB ./ GB(end);
    GA = GA ./ GA(end);
    
    GB = GB ./ GA(1);
    GA = GA ./ GA(1);

    if rtype == 0
        % Controllable canonical form
        F = zeros(length(GA)-1);
        F(end,:) = -GA(end:-1:2);
        F(1:end-1,2:end) = eye(length(GA)-2);
    
        L = zeros(length(GA)-1,1);
        L(end) = 1;
    
        H = zeros(1,length(GA)-1);
        H(1:length(GB)) = GB(end:-1:1);
    else
        % Observable canonical form
        F = zeros(length(GA)-1);
        F(:,end) = -GA(end:-1:2);
        F(2:end,1:end-1) = eye(length(GA)-2);
    
        L = zeros(length(GA)-1,1);
        L(1:length(GB)) = GB(end:-1:1);

        H = zeros(1,length(GA)-1);
        H(end) = 1;    
    end

    if nargout > 4
%        Pinf = lyap(F,F',L*q*L');
        Pinf = lyapchol(F,L*sqrt(q));
        Pinf = Pinf' * Pinf;
    end    
end


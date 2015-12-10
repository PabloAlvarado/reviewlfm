% [B,A] = se_taylor(n,s,ell)
%
% Taylor approximation of SE. If s or ell is symbolic,
% returns symbolic result.

% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

% Copyright (C) 2014 Simo Sarkka

function [B,A] = se_taylor(n,s,ell)

    if strcmp(class(s),'sym') || strcmp(class(ell),'sym')
        B = s^2 * sqrt(2*sym('pi')) * ell;
        A = sym(zeros(1,2*n+1));
    else
        B = s^2 * sqrt(2*pi) * ell;
        A = zeros(1,2*n+1);
    end
    
    i = 1;
    for k=n:-1:0
        A(i) = (ell^2/2)^k/factorial(k);
        i = i + 2;
    end
end


% [B,A] = se_power_pade(n,mpade,npade,s,ell)
%
% Pade approximant [mpade/npade]^n of SE.

% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

% Copyright (C) 2014 Simo Sarkka

function [B,A] = se_power_pade(n,mpade,npade,s,ell)

    [B,A] = se_pade(mpade,npade,s,ell);

    B = B ./ n.^((length(B)-1:-1:0)/2);
    A = A ./ n.^((length(A)-1:-1:0)/2);

    %S0 = polyval(B,0) / polyval(A,0);
    S0 = B(end) / A(end);
    B0 = B / S0;
    A0 = A;
    
    B = B0;
    A = A0;
    for k=1:n-1
        B = conv(B,B0);
        A = conv(A,A0);
%        B = my_conv(B,B0);
%        A = my_conv(A,A0);
    end
    
    B = S0 * B;
end


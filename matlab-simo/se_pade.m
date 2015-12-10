% [B,A] = se_pade(m,n,s,ell)
%
% Pade approximant [m/n] of SE. 

% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

% Copyright (C) 2014 Simo Sarkka

function [B,A] = se_pade(m,n,s,ell)
    order = m+n;
    N = order:-1:0;

    if strcmp(class(s),'sym')
        c = (ell^2/2).^N./factorial(N);
        c = sym(c);
    else
        c = (ell^2/2).^N./factorial(N);
%        c = exp(N .* log(ell^2/2) - gammaln(N+1));
    end

    [b,a] = pade_approx(c,n,m);

    A = zeros(1,2*length(a)-1);
    B = zeros(1,2*length(b)-1);
    
    if strcmp(class(s),'sym')
        A = sym(A);
        B = sym(B);
    end

    A(1:2:end) = a;
    B(1:2:end) = b;
 
    if strcmp(class(s),'sym')    
        B = B * s^2 * sqrt(2 * sym('pi')) * ell;
    else
        B = B * s^2 * sqrt(2 * pi) * ell;
    end
    
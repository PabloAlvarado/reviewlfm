%PADE_APPROX  Form Pade approximant for a series
%
% Syntax:
%   [a,b] = pade_approx(c,P,[Q,a0])
%
% Form a [P/Q] pade approximant from the given power series in c.
%
% Example:
%  order = 4;
%  N = order:-1:0;
%  P = 2;
%  Q = order - P;
%  c = (-1).^N./factorial(N);
%  [a,b] = pade_approx(c,P,Q);
%  x = 0:0.01:4;
%  plot(x,exp(-x),'-',x,polyval(c,x),':',x,polyval(b,x)./polyval(a,x),'--')
%  legend('Exact','Taylor','Pade')
%

% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

% Copyright (C) 2014 Simo Sarkka

function [a,b,G,y,P,Q] = pade_approx(c,P,Q,a0)

    if nargin < 3
        Q = [];
    end
    if nargin < 4
        a0 = [];
    end

    if isempty(Q)
        Q = length(c)-1-P;
    end
    if isempty(a0)
        a0 = 1;
    end

    
    if length(c) ~= P+Q+1
        error(sprintf('length(c) ~= P+Q+1 (c=%d,P=%d,Q=%d)',...
              length(c),P,Q));
    end    
    
    %
    % Solve a's
    %
    c = c(:)';
    G = zeros(Q,Q);
    if strcmp(class(c),'sym')
        G = sym(G);
    end
    y = -a0 * c(Q:-1:1)'; % -a0 [c{P+1},...,c{P+Q}]
    for i=1:Q
        ind1 = P+i-1;
        ind2 = P+i-Q;
        
        tmp = c((length(c)-ind1):(length(c)-max(ind2,0)));
        G(i,length(tmp):-1:1) = tmp;
    end
    a = [(G\y)' a0];
    
    %
    % Solve b's
    %
    b = conv(c,a);
%    b = my_conv(c,a);
    b = b(end-P:end);

    

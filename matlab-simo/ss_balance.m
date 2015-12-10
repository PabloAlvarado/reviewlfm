% [NF,NL,NH,T] = ss_balance(F,L,H)
%
% Given a state-space model
%
%  dx/dt = F x + L w
%      y = H x
%
% Find a matrix T such that the following equivalent
% state space model is numerically more stable than
% the original:
%
%  dz/dt = inv(T) F T z + inv(T) L w
%      y = H T z
%
% Returns NF = inv(T) F T, NL = inv(T) L, NH = H T.

% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

% Copyright (C) 2014 Simo Sarkka

function [NF,NL,NH,T] = ss_balance(F,L,H)
    [T,NF] = balance(F);
    NL = T\L;
    NH = H*T;
end


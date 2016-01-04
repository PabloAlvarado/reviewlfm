%
% Generate the sprint force and set parameters
%

  %
  % Parametrit
  %
  N = 5000;
  dt = 0.01;
  a = 0;
  b = 0.1;
  tt = (1:N)*dt;
  meas_step = 10;
  sd = 0.01;
  
  %
  % Unknown input
  %
  r1 = randn;
  r2 = randn;
  amp_ext = 0.5;
  u_ext = amp_ext*(sin(r1 + 0.43*tt) + sin(r1 + 0.19*tt));
  plot(tt,u_ext);
  
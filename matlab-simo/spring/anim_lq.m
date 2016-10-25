%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sprint with a basic LQ control
%
%                 u(t) = Kp (y(t) - a) + Kd dy/dt
%  d²y/dt² + b dy/dt + y = u(t) + a
%
%   dy_1/dt = y_2
%   dy_2/dt = -b y_2 - y_1 + a + Kp (y_1 - a) + Kd y_2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

  %%
  %
  rng(1);
  spring_param;

  %%
  % Animoi
  %
  clf;
  set(gcf,'DoubleBuffer','on');
  F_lq = [0 1; -1 -b];
  L_lq = [0;1];
  x = [-1;0];
  t = 0;
  X_lq = [];
  T_lq = [];
  G = lqr(F_lq,L_lq,diag([1 1]),0.05);
  Kp = -G(1);
  Kd = -G(2);
  for k=1:N
    u = Kp*(x(1) - a) + Kd*x(2);
    x = x + dt*(F_lq*x + L_lq*a + L_lq*u + L_lq*u_ext(k));
    t = t + dt;

    X_lq = [X_lq x];
    T_lq = [T_lq t];

    if rem(k,meas_step) == 0
        subplot(1,2,1);
        spring_plot(x,a);
        subplot(1,2,2);
        plot(T_lq,X_lq(1,:),[0 N*dt],[a a],'--');
        
        drawnow;
        pause(0.01);
    end
  end
  
  plot(T_lq,X_lq(1,:),[0 N*dt],[a a],'--');
  grid on;
  axis([0 N*dt -1 1]);

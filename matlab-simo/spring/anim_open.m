%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%
% Spring without a control
%
%   d²y/dt² + b dy/dt + y = a
%
%   dy_1/dt = y_2
%   dy_2/dt = -b y_2 - y_1 + a
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

  %%
  %
  rng(1);
  spring_param;

  %%
  % Animoi
  %
  if ~exist('quiet') || (quiet == 0)
      clf;
      set(gcf,'DoubleBuffer','on');
  end
  
  F = [0 1; -1 -b];
  L = [0;1];
  x = [-1;0];
  t = 0;
  X = [];
  Y = [];
  T = [];
  if ~exist('quiet') || (quiet == 0)
      clf;
  end
  
  for k=1:N
    x = x + dt*(F*x + L*a + L*u_ext(k));
    t = t + dt;
    
    if rem(k,meas_step) == 0
        y = x(1) + sd*randn;
        
        if ~exist('quiet') || (quiet == 0)
            subplot(1,2,1);
            spring_plot(x,a);
            subplot(1,2,2);
            plot(T,X(1,:),[0 N*dt],[a a],'--');

            drawnow;
            pause(0.01);
        end
    else
        y = NaN;
    end
    
    X = [X x];
    Y = [Y y];
    T = [T t];
  end
  
  ind = find(~isnan(Y));
  plot(T,X(1,:),T(ind),Y(ind),'.',[0 N*dt],[a a],'--');
  grid on;

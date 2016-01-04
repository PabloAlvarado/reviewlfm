function spring_plot(x,a)

  plot([0 0 2 2],[-3 3 3 -3],'b-',...
       [0.1 0.1 1.9 1.9 0.1],x(1)+[-0.2 0.2 0.2 -0.2 -0.2],'r-',...
       0.5*rem(1:10,2)+0.75,linspace(x(1)+0.2,3,10),'k-',...
       [0 2],[a a],'b--');
  axis off;
  axis([-2 4 -3.5 3.5]);
  % Define eigenfunction of the Laplacian in 1D
  % box with Dirichlet boundary conditions (x \in [-L, L])
  
  eigenval = @(n) (n(:)'*pi/2/L).^2;
  eigenfun = @(n,x) L^(-1/2)*sin(kron(n(:)'*pi,(x(:)+L)/2/L));

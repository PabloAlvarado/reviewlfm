function [eigenfun,eigenval,NN] = domain_cartesian(m,d,L)
%% domain_cartesian - Laplace operator eigendecomposition in a hypercube
% 
% Syntax:
%  [eigenfun,eigenval,NN] = domain_cartesian(m,d,L)
%
% In:
%   m  - Number of eigenfunctions
%   d  - Dimensionality
%   L  - Domain boundary [-L1,L1]x[-L2,L2]x...x[-Ln,Ln]
%      
% Out:
%   eigenfun - Function handle: eigenfun(n,x)
%   eigenval - Function handle: eigenval(n)
%   NN       - Indices to evaluate the handles at
% 
% Description:
%   This code returns the eigendecomposition of the Laplacian in
%   Cartesian coordinates (x1,x2,...) in [-L1,L1]x[-L2,L2]x...
%   with respect to indices (n1,n2,...). The function is vectorized
%   with respect to both the location x and the indices n.
%
% Copyright (c) 2014 Arno Solin
%
% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  % This is stupid, but at least we should get enough 
  % of basis function to choose from
  N = ceil(m^(1/d)*L/min(L));
  
  % Combined eigenfunction indices (checked the numbers)
  NN = ndgridm(N);

  % Define eigenvalues of the negative Laplacian in ND 
  % s.t. Dirichlet boundary. This forms an orthonormal basis.
  eigenval = @(n) sum((pi*bsxfun(@rdivide,n,2*L)).^2,2);
    
  % Sort and take only m most important eigenfunctions
  [~,ind] = sort(eigenval(NN)); NN = NN(ind(1:m),:);  

  % Define eigenfunction of the negative Laplacian in 2D 
  % s.t. Dirichlet boundary. This forms an orthonormal basis.
  eigenfun = @(n,x) laplace_eig_cart_dirichlet(n,x,L);

end

function [v]=laplace_eig_cart_dirichlet(n,x,L)
%% laplace_eig_cart_dirichlet - Laplace operator eigenfunctions in a hypercube
% 
% Syntax:
%  [v] = laplace_eig_cart_dirichlet(n,x,L)
%
% In:
%   n  - Eigenfunction indices
%   x  - Spatial locations [x1 x2]
%   L  - Domain boundary [-L1,L1]x[-L2,L2]x...x[-Ln,Ln]
%      
% Out:
%   v - The evaluated value
% 
% Description:
%   This code calculates the eigenvectors of the Laplacian in
%   Cartesian coordinates (x1,x2,...) in [-L1,L1]x[-L2,L2]x...
%   with respect to indices (n1,n2,...). The function is vectorized
%   with respect to both the location x and the indices n.
%
%   The corresponding eigenvalues can be calculated by
% 
%     eigenval = @(n) sum((pi*bsxfun(@rdivide,n,2*L)).^2,2);
%
% Copyright (C) 2012 Arno Solin
%

  % Allocate space
  v = zeros(size(x,1),size(n,1));

  % Evaluate eigenfunctions
  if size(x,2)==1
      for i=1:numel(n)
          v(:,i) = sqrt(1./L)*sin(pi*n(i)*(x(:)+L)/2/L);
      end
  else
      for i=1:size(n,1)
          % Eigenfunctions for x in Omega and n = (n1,n2,...nn), ni = 1,2,...,Ni
          v(:,i) = prod(bsxfun(@times,sqrt(1./L), ...
              sin(pi*bsxfun(@times,n(i,:)./L,bsxfun(@plus,x,L))/2)),2);
          if all(n(i,:)==0)
              v(:,i) = ones(size(x,1),1);
          end
      end
  end

end

function NN = ndgridm(N)
%% ndgridm - Expand index hypercude
%
% Syntax:
%  [NN] = ndgridm(N)
%
% In:
%   N  - Vector of max indices
%      
% Out:
%   NN - Matrix of index combinations
%
% Description:
%   A more felxible variant of 'ndgrid'. This functions gives combinations
%   of the indices in N such that, for example, for a 2D case we get
%   (1,1),(1,2),...,(1,N2),(2,1),...,(2,N2),(...),(N1,1),...,(N1,N2).
%   This function works for any dimension.
%
% Copyright (C) 2014 Arno Solin
%

  % Allocate space for indices
  NN = zeros(prod(N),numel(N));

  % For each level/diemsion
  if numel(N)==1

     % The lowest level
     NN(:,1) = (1:N)';
     
  else

    % This level
    n = 1:N(1);

    % Recursive call
    nn = ndgridm(N(2:end));

    % Assign values
    NN(:,1)     = kron(n,ones(1,prod(N(2:end))))';
    NN(:,2:end) = repmat(nn,[N(1) 1]);
   
  end

end

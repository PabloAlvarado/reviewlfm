function K = computeExpKernel(t, tp, ell)

N = length(t);
Np = length(tp);

K = zeros(N, Np);

for i=1:N
    for j=1:Np
        K(i,j) = exp(-ell*abs(t(i)- tp(j)));
    end
end
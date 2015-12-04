function K = computeLfmXExpKernel(t, tp, mass, damper, spring, sens, ell)

N = length(t);
Np = length(tp);

K = zeros(N, Np);

alpha = damper/(2*mass);
beta = sqrt((spring/mass) -  alpha^2);

for i=1:N
    for j=1:Np
        K(i,j) = computeLocalLfmXExpKernel(t(i), tp(j), ...
            alpha, beta, ell);
    end
end

C0 = sens/(mass*beta);
K = C0*K;

K = real(K);
end

function kl = computeLocalLfmXExpKernel(t, tp, alpha, beta, ell)

wl = (ell + alpha) - (1j)*beta;
w  = (ell + alpha) + (1j)*beta;
wt = (alpha - ell) - (1j)*beta;
wh = (alpha - ell) + (1j)*beta;


if t > tp
    kaux1 = exp(-alpha*t - ell*tp)*H2(0, t, wl, w, beta, tp);
    kaux2 = exp(-alpha*t + ell*tp)*H2(tp, t, wt, wh, beta, t);
    kl = (kaux1 + kaux2);
else
    kaux1 = exp(-alpha*t - ell*tp)*H2(0, t, wl, w, beta, t);
    kl = kaux1;
end

end

function vH2 = H2(tp, t, wt, wh, beta, tsup)
vH2 = (1/(2*(1j)*wt))*exp((1j)*beta*t)*exp(wt*tsup) ...
    - (1/(2*(1j)*wt))*exp((1j)*beta*t)*exp(wt*tp) ...
    - (1/(2*(1j)*wh))*exp(-(1j)*beta*t)*exp(wh*tsup) ...
    + (1/(2*(1j)*wh))*exp(-(1j)*beta*t)*exp(wh*tp);
end


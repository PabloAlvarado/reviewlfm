function K = computeLfmExpKernel(t, tp, mass, damper, spring, sens, ell)

N = length(t);
Np = length(tp);

K = zeros(N, Np);

alpha = damper/(2*mass);
beta = sqrt((spring/mass) -  alpha^2);

for i=1:N
    for j=1:Np
        K(i,j) = computeLocalLfmExpKernel(t(i), tp(j), ...
            alpha, beta, ell);
    end
end

C0 = sens^2/((mass^2)*(beta^2));
K = C0*K;

K = real(K);

end

function kl = computeLocalLfmExpKernel(t, tp, alpha, beta, ell)

wl = (ell + alpha) - (1j)*beta;
w  = (ell + alpha) + (1j)*beta;
wt = (alpha - ell) - (1j)*beta;
wh = (alpha - ell) + (1j)*beta;

if t > tp
    kaux1 = kpart(tp, t, wl, wt, wh, w, beta, tp);
    kaux2 = kpart(t, tp, wl, wt, wh, w, beta, tp);
    kaux3 = H2(tp, t, wt, wh, beta)*H2(0, tp, wl, w, beta);
else
    kaux1 = kpart(t, tp, wl, wt, wh, w, beta, t);
    kaux2 = kpart(tp, t, wl, wt, wh, w, beta, t);
    kaux3 = H2(t, tp, wt, wh, beta)*H2(0, t, wl, w, beta);
end

kl = exp(-alpha*(t + tp))*(kaux1 + kaux2 + kaux3);

end

function kregion = kpart(tp, t, wl, wt, wh, w, beta, tsup) 
kregion = - H0(tp, t, wl, wt, beta, beta, tsup) + H0(tp, t, wl, wh, beta, -beta, tsup) ...
    + H1(tp, t, wl, wt, beta, beta, tsup) - H1(tp, t, wl, wh, beta, -beta, tsup) ...
    + H0(tp, t, w, wt, -beta, beta, tsup) - H0(tp, t, w, wh, -beta, -beta, tsup) ...
    - H1(tp, t, w, wt, -beta, beta, tsup) + H1(tp, t, w, wh, -beta, -beta, tsup);

end

function  vH0 = H0(tp, t, wl, wt, bdp, bd, tsup)
vH0 = (1/(4*wl))*(1/(wt + wl))*exp((1j)*bdp*tp)*exp((1j)*bd*t)*...
    (exp((wt + wl)*tsup)-1);
end

function vH1 = H1(tp, t, wl, wt, bdp, bd, tsup)
vH1 = (1/(4*wl*wt))*exp((1j)*bdp*tp)*exp((1j)*bd*t)*(exp(wt*tsup)-1);
end

function vH2 = H2(tp, t, wt, wh, beta)
vH2 = (1/(2*(1j)*wt))*exp((1j)*beta*t)*exp(wt*t) ...
    - (1/(2*(1j)*wt))*exp((1j)*beta*t)*exp(wt*tp) ...
    - (1/(2*(1j)*wh))*exp(-(1j)*beta*t)*exp(wh*t) ...
    + (1/(2*(1j)*wh))*exp(-(1j)*beta*t)*exp(wh*tp);
end




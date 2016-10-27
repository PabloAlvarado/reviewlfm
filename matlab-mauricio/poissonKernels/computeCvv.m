function const = computeCvv(sigmax, lengthX, gamma, pwz1, pwz2, n, m)

% COMPUTECVV Compute a constant needed for the Poisson kernel
%
% COPYRIGHT : Mauricio A. Alvarez, 2016

% KERN

if n == m    
    Wox = exp((gamma(n)*sigmax/2)^2)*(Faddeeva_erfc(pwz1(n)) - Faddeeva_erfc(pwz2(n)));
    const = (sigmax*sqrt(pi)*lengthX/2)*(real(Wox) ...
        - imag(Wox)*((sigmax^2*n*pi)/(2*lengthX^2) + (1/(n*pi)))) ...
        +(sigmax^2/2)*(exp(-(lengthX/sigmax)^2)*cos(n*pi) - 1);
else
    if mod(n+m,2)==1
        const = 0;
    else
        Woxm = exp((gamma(m)*sigmax/2)^2)*(Faddeeva_erfc(pwz1(m)) - Faddeeva_erfc(pwz2(m)));
        Woxn = exp((gamma(n)*sigmax/2)^2)*(Faddeeva_erfc(pwz1(n)) - Faddeeva_erfc(pwz2(n)));
        const = ((sigmax*lengthX)/(sqrt(pi)*(m^2-n^2)))*(n*imag(Woxm) - m*imag(Woxn));
    end
end

function const = computeCvv(sigmax, lengthX, gamma, pwz1, pwz2, n, m)

% COMPUTECVV Compute a constant needed for the Poisson kernel
%
% COPYRIGHT : Mauricio A. Alvarez, 2016

% KERN

if n == m
    Wox = pwz1(n) - exp(-(lengthX/sigmax)^2)*exp(-gamma(n)*lengthX)*pwz2(n);
    const = (sigmax*sqrt(pi)*lengthX/2)*(real(Wox) ...
        - imag(Wox)*((sigmax^2*n*pi)/(2*lengthX^2) + (1/(n*pi)))) ...
        +(sigmax^2/2)*(exp(-(lengthX/sigmax)^2)*cos(n*pi) - 1);
else
    if mod(n+m,2)==1
        const = 0;
    else
        Woxm = pwz1(m) - exp(-(lengthX/sigmax)^2)*exp(-gamma(m)*lengthX)*pwz2(m);
        Woxn = pwz1(n) - exp(-(lengthX/sigmax)^2)*exp(-gamma(n)*lengthX)*pwz2(n);
        const = ((sigmax*lengthX)/(sqrt(pi)*(m^2-n^2)))*(n*imag(Woxm) - m*imag(Woxn));
    end
end

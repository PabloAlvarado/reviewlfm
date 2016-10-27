clc
close all
clear

x = 0.1;
xp = 0.5;
y = 0.6;
yp = 0.9;
nterms = 5;
sigmax = 0.3;
sigmax2 = sigmax^2;
sigmay = 0.3;
sigmay2 = sigmay^2;
a = 2;
b = 2;

C = 16/((a*b)^2);

knum = 0;
for n=1:nterms
    for np = 1:nterms
        for m = 1:nterms
            for mp = 1:nterms
                pn = (n*pi)/a;
                pnp = (np*pi)/a;
                qm = (m*pi)/b;
                qmp = (mp*pi)/b;
                pnqm = pn^2 + qm^2;
                pnpqmp = pnp^2 + qmp^2;
                kterm = (1/(pnqm*pnpqmp))*sin(pn*x)*sin(pnp*xp)*...
                    sin(qm*y)*sin(qmp*yp)*...
                    integral2(@(xi,xip) sin(pn*xi).*sin(pnp*xip).* ...
                    exp(-((xi -xip).^2)./sigmax2),0,a,0,a, 'RelTol', 1e-10)*...
                    integral2(@(eta,etap) sin(qm*eta).*sin(qmp*etap).* ...
                    exp(-((eta-etap).^2)./sigmay2),0,b,0,b, 'RelTol', 1e-10);
                knum = knum + kterm;
            end
        end
    end
end
knum = C*knum;
%%% Now using the kernel
kern = kernCreate([x y], 'poisson');
kern.nTerms = nterms;
kern.inverseWidthX = 2/sigmax2;
kern.inverseWidthY = 2/sigmay2;
kern.lengthX = a;
kern.lengthY = b;
kanal = poissonKernCompute(kern, [x y], [xp yp]);









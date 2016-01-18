clc
clear
close all
s = RandStream('mt19937ar', 'Seed', 1e4);
RandStream.setGlobalStream(s);

n = 4;
nX = 50;
nY = 50;
xi = linspace(-1,1, nX)';
yi = linspace(-1,1, nY)';
[xx1, xx2] = meshgrid(xi, yi);
X = [xx1(:) xx2(:)];

lengthX = max(xi) - min(xi);
lengthY = max(yi) - min(yi);

aij = 100*randn(n, n);
bi = 10*rand(n, 1);
cj = 10*rand(n, 1);
pi = 5*rand(n, 1);
qj = 5*rand(n, 1);

isexp = false;

ff_p = zeros(nX*nY, 1);
uu = zeros(nX*nY, 1);

for i=1:n
    for j=1:n
        if isexp            
            ff_p = ff_p  - (aij(i, j)/(bi(i)^2 + cj(j)^2))*...
                exp(bi(i)*X(:,1) + cj(j)*X(:, 2));
        else            
            ff_p = ff_p + (aij(i, j)/(bi(i)^2 + cj(j)^2))*...
                sin(bi(i)*X(:,1) + pi(i)).*sin(cj(j)*X(:, 2) + qj(j));
        end
    end
end

pcolor(xx1, xx2, reshape(ff_p, size(xx1)))
colorbar
figure
surf(xx1, xx2, reshape(ff_p, size(xx1)))

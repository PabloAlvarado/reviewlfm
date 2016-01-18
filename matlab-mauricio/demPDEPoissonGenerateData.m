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
            uu = uu + aij(i, j)*exp(bi(i)*X(:,1) + cj(j)*X(:, 2));
            ff_p = ff_p  - (aij(i, j)/(bi(i)^2 + cj(j)^2))*...
                exp(bi(i)*X(:,1) + cj(j)*X(:, 2));
        else
            uu = uu + aij(i, j)*sin(bi(i)*X(:,1) + pi(i))...
                .*sin(cj(j)*X(:, 2) + qj(j));
            ff_p = ff_p + (aij(i, j)/(bi(i)^2 + cj(j)^2))*...
                sin(bi(i)*X(:,1) + pi(i)).*sin(cj(j)*X(:, 2) + qj(j));
        end
    end
end
figure
pcolor(xx1, xx2, reshape(ff_p, size(xx1)))
colorbar
figure
surf(xx1, xx2, reshape(ff_p, size(xx1)))
figure
pcolor(xx1, xx2, reshape(uu, size(xx1)))
colorbar
figure
surf(xx1, xx2, reshape(uu, size(xx1)))
savedata = false;
if savedata
    save('./dataPoissonExample2.mat', 'xx1', 'xx2', 'uu', 'ff_p')
end

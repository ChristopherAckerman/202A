beta = 0.99;
epsilon = 1e-6;
alpha = 0.36;
delta = .02;

ks = ((1 - beta*(1 - delta))/(alpha * beta))^(1/(alpha  - 1));

kmin = 0;
kmax = 50;
grid_size = 10000;

dk = (kmax-kmin)/(grid_size - 1);
kgrid = linspace(kmin, kmax, grid_size);
v = zeros(grid_size, 1);
dr = zeros(grid_size, 1);
norm = 1;

while norm > epsilon;
    for i=1:grid_size
        tmp = (kgrid(i)^alpha + (1 - delta)*kgrid(i) - kmin);
        imax = min(floor(tmp/dk) + 1, grid_size);
        
        c = kgrid(i)^alpha + (1 - delta)*kgrid(i) - kgrid(1:imax);
        util = log(c);
        [tv(i), dr(i)] = max(util + beta*(1:imax));
    end;
    norm = max(abs(tv - v));
    v = tv;
end;
clear;
m = 2; % number of pre-determined variables: k, z
n = 1; % number of forward-looking variables: c

% values for calibration
barR = 1.05;
ItoYratio = 0.2;
RKtoYratio = 0.3;
CtoYratio = 1 - ItoYratio;

% calibration
bbeta = 1 / barR;
ttheta = RKtoYratio;
ggamma = 0.02;
ddelta = 1 - ((1 + ggamma) * (CtoYratio - 1 + bbeta * ttheta)) / (bbeta * ttheta + CtoYratio * bbeta - bbeta);

% steady state
Kss = (bbeta * ttheta / ((1 + ggamma) - bbeta * (1 - ddelta)))^(1 - ttheta);
Yss = Kss^ttheta;
Css = Yss + (1 - ddelta - (1 + ggamma)) * Kss;
Rss = barR;

% define matrices

rrho = 0.95;
M31 = 0;
M32 = 0;
M33 = rrho;

M23 = Kss^(ttheta - 1) / (1 + ggamma);
M22 = (ttheta * Kss^(ttheta - 1) + 1 - ddelta) / (1 + ggamma);
M21 = -Css / ((1 + ggamma) * Kss);

tmp = (bbeta * ttheta * Kss^(ttheta - 1) / (1 + ggamma));
M13 = tmp * ((ttheta - 1) * M23 + M33);
M12 = tmp * (ttheta - 1) * (ttheta * Kss^(ttheta - 1) + 1 - ddelta) / (1 + ggamma);
M11 = 1 - tmp * (ttheta - 1) * Css / ((1 + ggamma) * Kss);

% diagonalize
M = [M11, M12, M13; M21, M22, M23; M31, M32, M33];
[Gamma, Lambda] = eig(M);
Lambda = diag(Lambda);
% sort eigenvalues in ascending order
[unused, order] = sort(abs(Lambda), 'ascend');
% reorder J and make diagonal again
Lambda = diag(Lambda(order));
% reorder eigenvectors
Gamma = Gamma(:, order);

% check number of eigenvalues outside unit circle equal to n
if(sum(abs(diag(Lambda)) > 1) ~= n)
    return;
end

if(sum(abs(diag(Lambda)) > 1) == n)
    disp("By Blanchard Khan proposition 1, we have a unique solution.");
end

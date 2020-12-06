clear;
m = 2; % number of pre-determined variables: k, z
n = 1; % number of forward-looking variables: c

% values for calibration
barR = 1.05;
ItoYratio = 0.2;
RKtoYratio = 0.3;
ybar = 1;
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

% partition matrices
Gammainv = inv(Gamma);
G11 = Gammainv(1 : m, 1 : n);
G12 = Gammainv(1 : m, n + 1 : m + n);
G21 = Gammainv(m + 1 : m + n, 1 : n);
G22 = Gammainv(m + 1 : m + n, n + 1 : m + n);
Lambda1 = Lambda(1 : m , 1 : m);
Lambda2 = Lambda(m + 1 : m + n, m + 1 : m + n);

% state variables solution
% E x_t+t = H*x_t
H = inv(-G11 * inv(G21) * G22 + G12) * Lambda1 * (-G11 * inv(G21) * G22+G12);
% Jump Variables : Policy function
THETA = -inv(G21) * G22;

% Simulation
Mtilde11 = M11;
Mtilde12 =[M12, M13];
Mtilde21 = [M21 ;M31];
Mtilde22 = [M22, M23; M32, M33];
A12=Mtilde11 * THETA + Mtilde12;
A22=Mtilde21 * THETA + Mtilde22;
A=[zeros(n, size(n, n)), A12 ; 
   zeros(m, n), A22 ];
T =100; % number of periods
x2v(1, 1) = 0;
x2v(1, 2) = 0;
x1v(1, 1) = THETA * x2v (1 ,:)';
X = zeros(T, m + n);
error = normrnd(0, 0.01, [T, 1]);
X(1, :) = [x1v ,x2v];
for i = 2 : 1 : T
X(i ,:) = (A * X(i - 1, :)' + [0, 0, error(i)]')';
end
GDP = Yss * exp((X(:, 2) * ttheta + X(:, 3)));
plot(linspace(0, T, T), GDP, 'LineWidth', 3, 'Color', 'k')
grid on
xlabel('Time')
ylabel('GDP t')
hold on
plot(linspace(0, T, T), GDP * 0 + Yss, 'LineWidth', 3, 'Color', 'r')
legend('GDP (levels)', 'Yss')

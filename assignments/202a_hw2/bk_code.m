clear;
% Neoclassical Model 
% Deterministic Model:
m =1; % number of pre - determined variables
n =1; % number of forward - looking variables
% values for calibration
barR       = 1.07;
ItoYratio  = 21/100;
RKtoYratio = 1-60/100;
ybar       = 1;
% Calibration
ssigma  = 1     ;
bbeta   = 1/barR;
ttheta  = RKtoYratio;
ddelta  = (1-barR)/(1-ttheta/ItoYratio);
Abar    = ((barR-1+ddelta)/ttheta)^(ttheta)*ybar^((1-ttheta));
% Steady State
Kss     = (Abar*ttheta/(barR-1+ddelta))^(1/(1-ttheta));
Yss     =  Abar*Kss^ttheta;
Iss     =  ddelta*Kss;
Css     =  Yss -Iss;
Rss     =  barR;
% Define matrices
% c_t+1 = [ M11    M12    ] [ c_t k_t]
% k_t+1 = [ M21    M22    ] [ c_t k_t]
M22= Abar*ttheta*Kss^(ttheta-1)+(1-ddelta); 
M21= -Css/Kss;
M12= bbeta*ttheta*Abar*Kss^(ttheta-1)*((ttheta-1)*M22);
M11= 1+bbeta*ttheta*Abar*Kss^(ttheta-1)*(ttheta-1)*M21;
% Diagonalize
M=[M11,M12; M21,M22];
[Gamma , Lambda ]= eig ( M );
Lambda = diag ( Lambda );
[ unused, order ]= sort ( abs ( Lambda) ,'ascend'); % sort eigenvalues in ascending order
Lambda = diag ( Lambda ( order )); % reorder J and make diagonal again
Gamma= Gamma(: , order ); % reorder eigenvectors

% check number of eigenvalues outside unit circle equal to n
if ( sum ( abs ( diag ( Lambda )) >1)~= n )
 return;
end

% partition matrices
Gammainv = inv ( Gamma );
G11 = Gammainv (1: m ,1: n );
G12 = Gammainv (1: m , n +1: m + n );
G21 = Gammainv ( m+1: m +n ,1: n );
G22 = Gammainv ( m+1: m +n , n +1: m + n );
Lambda1 = Lambda (1: m ,1: m );
Lambda2 = Lambda ( m +1: m +n , m +1: m + n );
% State Variables solution
% E x_t+1 = H*x_t
H=inv(-G11*inv(G21)*G22+G12)*Lambda1*(-G11*inv(G21)*G22+G12);
% Jump Variables: Policy function
THETA=-inv(G21)*G22;
% Simulation
T =100; % number of periods
x1 = zeros (n , T);
x2 = zeros (m , T );
x2(1,1)=0.2;
for i =1:1:T
x2 (: , i +1)= H * x2 (: , i );
x1 (: , i ) =  THETA * x2 (: , i );
end
% plot evolution of capital
plot(linspace(0,T,T+1),Kss*exp(x2),'LineWidth',3,'Color','k')
grid on
xlabel('Time')
ylabel('K_t')
hold on
plot(linspace(0,T,T+1),Kss*exp(x2)*0+Kss,'LineWidth',3,'Color','r')
legend('Capital','Kss')

% plot evolution of consumption
plot(linspace(0,T,T),Css*exp(x1),'LineWidth',3,'Color','k')
grid on
xlabel('Time')
ylabel('K_t')
hold on
plot(linspace(0,T,T),Css*exp(x1)*0+Css,'LineWidth',3,'Color','r')
legend('Consumption','Css')

%% Avoiding numerical error
A12=M11*THETA+M12;
A22=M21*THETA+M22;
A=[zeros(n,size(A12,2)),A12;zeros(m,size(A22,2)),A22];
T =100; % number of periods
x2v(1,1)=0.2;
x1v(1,1)=THETA*x2v(:,1);
X=zeros(T,m+n);
X(1,:)=[x1v,x2v];
for i =2:1:T
X(i,:)=(A*X(i-1,:)')';
end
% plot evolution of capital
plot(linspace(0,T,T),Kss*exp(X(:,2)),'LineWidth',3,'Color','k')
grid on
xlabel('Time')
ylabel('K_t')
hold on
plot(linspace(0,T,T),Kss*exp(X(:,2))*0+Kss,'LineWidth',3,'Color','r')
legend('Capital','Kss')

% plot evolution of consumption
plot(linspace(0,T,T),Css*exp(X(:,1)),'LineWidth',3,'Color','k')
grid on
xlabel('Time')
ylabel('C_t')
hold on
plot(linspace(0,T,T),Css*exp(x1)*0+Css,'LineWidth',3,'Color','r')
legend('Consumption','Css')

%% Stochastic Model
rrho=0.9;
M31=0;
M32=0;
M33=rrho;
M23= Abar*Kss^(ttheta-1) ;
M13= bbeta*Abar*Kss^(ttheta-1)*(M33+(ttheta-1)*M23);
M=[M11,M12,M13; M21, M22, M23; 0,0,M33];
m =2; % number of pre - determined variables
n =1; % number of forward - looking variables
[Gamma , Lambda ]= eig ( M );
Lambda = diag ( Lambda );
[ unused, order ]= sort ( abs ( Lambda) ,'ascend'); % sort eigenvalues in ascending order
Lambda = diag ( Lambda ( order )); % reorder J and make diagonal again
Gamma= Gamma(: , order ); % reorder eigenvectors

% check number of eigenvalues outside unit circle equal to n
if ( sum ( abs ( diag ( Lambda )) >1)~= n )
 return;
end

% partition matrices
Gammainv = inv ( Gamma );
G11 = Gammainv (1: m ,1: n );
G12 = Gammainv (1: m , n +1: m + n );
G21 = Gammainv ( m+1: m +n ,1: n );
G22 = Gammainv ( m+1: m +n , n +1: m + n );
Lambda1 = Lambda (1: m ,1: m );
Lambda2 = Lambda ( m +1: m +n , m +1: m + n );
% State Variables solution
% E x_t+1 = H*x_t
H=inv(-G11*inv(G21)*G22+G12)*Lambda1*(-G11*inv(G21)*G22+G12);
% Jump Variables: Policy function
THETA=-inv(G21)*G22;
% Simulation
T =20; % number of periods
Mtilde11=M11;
Mtilde12=[M12,M13];
Mtilde21=[M21;M31];
Mtilde22=[M22,M23;M32,M33];
A12=Mtilde11*THETA+Mtilde12;
A22=Mtilde21*THETA+Mtilde22;
A=[zeros(n,size(n,n)),A12;zeros(m,n),A22];
T =100; % number of periods
x2v(1,1)=0;
x2v(1,2)=0;
x1v(1,1)=THETA*x2v(1,:)';
X=zeros(T,m+n);
error=normrnd(0,0.01,[T,1]);
X(1,:)=[x1v,x2v];
for i =2:1:T
X(i,:)=(A*X(i-1,:)'+[0,0,error(i)]')';
end
GDP=Yss*exp((X(:,2)*ttheta+X(:,3)));
plot(linspace(0,T,T),GDP,'LineWidth',3,'Color','k')
grid on
xlabel('Time')
ylabel('GDP_t')
hold on
plot(linspace(0,T,T),GDP*0+Yss,'LineWidth',3,'Color','r')
legend('GDP (levels)','Yss')

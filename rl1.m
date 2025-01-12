clear

n=4;
m=2;
A=[2 1 3 0; 2 -1 0 1; 1 -2 -2 1; 0 1 1 -1];
%eig(A)

B=1*[ 1 0 ; 2 -1; 3 2; 0 1];
    
C=[ 1 3 -1 0; 2 -1 0 1];

D=[ 1 2; 0 0];

sys=ss(A,B,C,D);
z=tzero(sys);
A=A-z(3)*eye(n);
sys=ss(A,B,C,D);
eig(A)
z=tzero(sys)

j=complex(0,1);
A0=A-j*eye(n);
D0=D-C*A0^(-1)*B;
eig(D0);

sys0=ss(A0^(-1),A0^(-1)*B, C*A0^(-1), -D0)
z0=tzero(sys0);

[X,K,L] = icare(A0^(-1),A0^(-1)*B,(A0^(-1))'*C'*C*A0^(-1),D0'*D0,-(A0^(-1))'*C'*D0)
X=real(X);

% proverka
MLLM=[A'*X+X*A+C'*C, X*B+C'*D; B'*X+D'*C, D'*D]
e=eig(MLLM)
% %%%%%%%%%%%%%%%%%%

[VV,DD]=eig(MLLM);
DD1=sqrtm(DD(n+1:n+m,n+1:n+m));
ML=VV(:,n+1:n+m)*DD1;
ML=ML';
M=ML(:,1:n);
L=ML(:,n+1:n+m);

% proverka
MLLM-[M'; L']*[M,L]
% %%%%%%%%%%%%%%%%%%%

% proverka
EE=[eye(n) zeros(n,m); zeros(m,n+m)];
AA=[A B ; M L];
ee=eig(AA,EE);
% %%%%%%

pom2=null3(L');  % null3 e moja
pom1=null(pom2');
Up=[pom1,pom2];
bL=pom1'*L;
bM=pom1'*M;
tM=pom2'*M;

% proverka
pom=[bM, bL; tM, zeros(1,2)];
pom'*pom-[M'; L']*[M,L]
% %%%%%%%%%%%%%%%%%%%%%%%%%%

% proverka
EE=[eye(n) zeros(n,m); zeros(m,n+m)];
AA=[A B ; bM bL; tM 0 0];
eee=eig(AA,EE);
% %%%%%%

[UU,SS,VV]=svd(bL);
V=VV*[(SS(1,1))^(-1) 0; 0 1];
bL=bL*V;
BB=B*V;
B1=BB(:,1);
B2=BB(:,2);
D1=D*V(:,1);
D2=D*V(:,2);
A1=A-B1*bM;

% proverka
EE=[eye(n) zeros(n,m); zeros(m,n+m)];
AA=[A B1 B2; bM bL; tM 0 0];
eeee=eig(AA,EE);
% %%%%%%

% simulacija

x0=[1; 2; -1; 3]; % poceten uslov
Jmin=x0'*X*x0;

N=100000;
dt=0.000000001;

for i=1:N
t(i)=i*dt;
end

B20=null(B2');
cE=[B20'; zeros(1,n)];
cF=[B20';  tM];
Skok=cF^(-1)*cE;

x0pl=Skok*x0; 
tM0=null(tM);
%Aa=(B20'*tM0)^(-1)*B20'*A1*tM0;
Aa=A-B1*bM-B2*(tM*B2)^(-1)*tM*A1;
%pom=pinv(tM0)*x0pl;
pom=x0pl;
pomy=(C-D1*bM)*x0pl;

uu1(1)=-bM*x0pl; 
uu2(1)=-(tM*B2)^(-1)*tM*A1*x0pl;
xx(:,1)=x0pl;
for i=1:N-1

    pom=expm(Aa*i*dt)*pom;
    %xx(:,i+1)=tM0*pom;
    xx(:,i+1)=pom;
    yy(:,i+1)=(C-D1*bM)*xx(:,i+1);
    
    uu1(i+1)=-bM*xx(:,i+1);
    uu2(i+1)=-(tM*B2)^(-1)*tM*A1*xx(:,i+1);
  
end

plot(t,yy(1,:),t,yy(2,:)) % plot for Figure 2

plot(t,uu1,t,uu2) % plot for Figure 3

% suboptimal control
al=-200;
Ab=[B20'; tM]^(-1)*[B20'*A1; al*tM];
Jms=0;
h(1)=tM*x0;
xxx(:,1)=x0;
uuu(1,1)=-bM*x0;
uuu(2,1)=(tM*B2)^(-1)*(-tM*A1*x0+al*h(1));
for i=1:N-1
    xxx(:,i+1)=expm(Ab*i*dt)*xxx(:,i);
    h(i+1)=tM*xxx(:,i+1);
    uuu(1,i+1)=-bM*xxx(:,i+1);
    
    
    uuu(2,i+1)=(tM*B2)^(-1)*(-tM*A1*xxx(:,i+1)+al*h(i+1));
    yyy(:,i+1)=(C-D1*bM)*xxx(:,i+1);
end

plot(t,uuu(1,:),t,uuu(2,:)) % plot for Figure 5

Cb=(C-D1*bM);
R2=lyap(Ab',Cb'*Cb);
x0'*R2*x0 % ova e priblizhno isto so Jmin


dal=0.01;
for i=1:20000
    all(i)=-(i+500)*dal;
    Abl=[B20'; tM]^(-1)*[B20'*A1; all(i)*tM];
    R2l=lyap(Abl',Cb'*Cb);
    crit(i)=x0'*R2l*x0;
end

plot(-all,crit) % Plot for Figure 4


% Comparison of methods for solving the LMI
j=complex(0,1);
A0=(A-j*eye(n))^(-1);
D0=D-C*A0*B;
eig(D0);
B0=-A0*B;
C0=C*A0;

sys0=ss(A0,B0, C0, D0);
z0=tzero(sys0);

% The first method is using icare. Output XXr
[XX,KK,LL] = icare(A0,B0,C0'*C0,D0'*D0,C0'*D0);
XXr=real(XX);

% The first method is the eigenvector method. Output Yr
AA=A0-B0*(D0'*D0)^(-1)*D0'*C0;
QQ=C0'*C0-C0'*D0*(D0'*D0)^(-1)*D0'*C0;
SS=B0*(D0'*D0)^(-1)*B0';
H=[AA -SS; -QQ -AA'];
e=eig(H);

v1=null(H-e(4)*eye(2*n));
v2=null(H-e(6)*eye(2*n));
v3=null(H-e(7)*eye(2*n));
v4=null(H-e(8)*eye(2*n));
v=[v1,v2,v3,v4];
U=v(1:n,:);
V=v(n+1:2*n,:);
Y=V*U^(-1);
Yr=real(Y);

% plotting xx needed for Figure 1
numPoints = 50; % Number of points to be sampled
indices = round(linspace(1, length(t), numPoints));

% Sampled data
t_sampled = t(indices);
xx_sampled = xx(:, indices);
plot(t_sampled, xx_sampled(1,:), 'o', 'LineStyle', 'none', 'MarkerSize', 2, 'Color', 'k');
hold on;

% Add the remaining plots
plot(t, xx(2,:), '-','Color','b');
plot(t, xx(3,:), '--','Color','r');
plot(t, xx(4,:), '-.','Color','g','LineWidth',1);

% Hold off to finish the plot
hold off;

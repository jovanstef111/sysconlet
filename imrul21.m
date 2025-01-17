clear
A=[3 0 -2 2 0;
1 -3 2 -1 5;
-2 8 3 -1 -8;
-5 3 2 -2 -4
1 -5 0 0 6];
B1=[1 0;
0 1;
-1 -1;
-2 -1;
0 1];
B2=[-1 0;
0 1;
0 -2;
1 -1;
-1 2];
Q=[18 -4 0 9 13;
-4 15 8 -6 -5;
0 8 6 -3 1;
9 -6 -3 6 6;
13 -5 1 6 13];
S2=[-3 -6;
9 2;
3 2;
-3 -4;
-6 -2];
hR=[9 0;
0 4];

%%% C,D
C=[-2 1 0 -1 -2;
-2 -2 -2 0 -2;
-1 3 1 -1 -2;
-3 1 1 -2 -1];
D1=zeros(4,2);
D2=[0 0;
0 0;
3 0;
0 2];

% solution in the Imrul's paper
Kmax=9.6*[ 4  2  0  2  0;
    2 1 0 1 0;
    0 0 0 0 0;
    2 1 0 1 0;
    0 0 0 0 0]

% invariant subspace method
B=[B1,B2];
D=[D1,D2];
n=5;

j=complex(0,1);
A0=(A+j*eye(n))^(-1);
D0=D-C*A0*B;
eig(D0);
B0=-A0*B;
C0=C*A0;


AA=A0-B0*(D0'*D0)^(-1)*D0'*C0;
SS=B0*(D0'*D0)^(-1)*B0';

[VV,DD]=eig(-AA');
UU=lyap(AA,-DD,-SS*VV);
XXX=VV*UU^(-1);
XXXr=real(XXX)

Kmax-XXXr


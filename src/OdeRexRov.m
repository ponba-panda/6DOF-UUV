%Importing the Model Matrices
clc;
clear all;
ModelRexRov
Xinitial = [0 0 0 0 0 0 10 0 0 0 0 0].';
%options = odeset('RelTol',1e-3,'AbsTol',1e-2);
Jm = matlabFunction(J,'File','Jmfile','Vars',{argnames(J).'});
Jinv = inv(J);
Jminv = matlabFunction(Jinv,'File','Jminvfile','Vars',{argnames(Jinv).'});
gm = matlabFunction(g,'File','gmfile','Vars',{argnames(g).'});
Dm = matlabFunction(D,'File','Dmfile','Vars',{argnames(D).'});
Cm = matlabFunction(C,'File','Cmfile','Vars',{argnames(C).'});
Minv = inv(M);
%[t Xfinal] = ode15s(@(t,State) odefun(t,State,M,Cm,Jm,Dm,gm),[0 100],Xinitial) 

%ODE Function
function Xdot = odefun(t,State,M,Cm,Jm,Dm,gm)
%global M C J D g;
%Define state X = [eta V]
%X = [eta;V]; %12 dimensions
x = State(1);
y = State(2);
z = State(3);
phi = State(4);
theta = State(5);
psi = State(6);
u = State(7);
v = State(8);
w = State(9);
p = State(10);
q = State(11);
r = State(12);
Xdot = zeros(12,1);
%M1 is substituted M
%M1 = M;
X1 = [x,y,z,phi,theta,psi,u,v,w,p,q,r].';
J1 = vpa(Jm(X1(1:6)),5);
%g1 = vpa(gm(X1(1:6)),5);
g1 = zeros(6,1);
D1 = vpa(Dm(X1(7:12)),5);
C1 = vpa(Cm(X1(7:12)),5);
Xdot(1:6) = J1*X1(7:12);
Xdot(7:12) = M\(-g1 -D1*X1(7:12)-C1*X1(7:12));
end
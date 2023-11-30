clc;
clear all;
syms V v1 v2 u v w p q r x y z theta phi psi eta;
%global M C D g J;
syms C(u,v,w,p,q,r) D(u,v,w,p,q,r) g(x,y,z,phi,theta,psi) J(x,y,z,phi,theta,psi);
v1 = [u v w].';
v2 = [p q r].';
V = [v1.' v2.'].';
eta = [x y z phi theta psi].';

B = 18484.168; %Buoyant Force
W = 18274.75; %Weight
% Body frame origin coinciding with cog => rbg = [0 0 0].'
S = @skew;
xb = 0.0822;
yb = -0.0077;
zb = 0.3872;
rb = [xb yb zb].'; %Centre of Buoyancy
MRB = [
1862.87 0 0 0 0 0;
0 1862.87 0 0 0 0;
0 0 1862.87 0 0 0;
0 0 0 525.39 1.44 33.41;
0 0 0 1.44 794.20 2.60;
0 0 0 33.41 2.60 691.23;
];
MRB_dec = mat2cell(MRB,[3 3],[3 3] );

MA = [
779.79 -6.8773 -103.32 8.5426 -165.54 -7.8033;
-6.8773 1222 51.29 409.44 -5.8488 62.726;
-103.32 51.29 3659.9 6.1112 -386.42 10.774;
8.5426 409.44 6.1112 534.9 -10.027 21.019;
-165.54 -5.8488 -386.42 -10.027 842.69 -1.1162;
-7.8033 62.726 10.775 21.019 -1.1162 224.32;
];
MA_dec = mat2cell(MA,[3 3],[3 3] );

CRB = [
MRB_dec{1,1}*S(v2) S(v2)*MRB_dec{1,2};
MRB_dec{2,1}*S(v2) -S(MRB_dec{2,2}*v2);
];
CRB = vpa(CRB,5);

CA = [
zeros(3) -S(MA_dec{1,1}*v1 + MA_dec{1,2}*v2);
-S(MA_dec{1,1}*v1+MA_dec{1,2}*v2) -S(MA_dec{2,1}*v1+MA_dec{2,2}*v2);
];
CA = vpa(CA,5);
DL = vpa(diag([74.82 69.48 728.40 268.80 309.77 105.00]),10);
DNL = vpa(diag([748.22*abs(u) 992.53*abs(v) 1821.01*abs(w) 672*abs(p) 774.44*abs(q) 523.27*abs(r)]),10);
D(u,v,w,p,q,r) = vpa(DL+DNL,5);

g(x,y,z,phi,theta,psi) = vpa([
(W-B)*sin(theta);
-(W-B)*cos(theta)*sin(phi);
-(W-B)*cos(theta)*cos(phi);
-(-yb*B)*cos(theta)*cos(phi)+(-zb*B)*cos(theta)*sin(phi);
(-zb*B)*sin(theta)+(-xb*B)*cos(theta)*cos(phi);
-(-xb*B)*cos(theta)*sin(phi) - (-yb*B)*sin(theta);
],5);

Rnb = [
cos(psi)*cos(theta) -sin(psi)*cos(phi)+cos(psi)*sin(theta)*sin(phi) sin(psi)*sin(phi)+cos(psi)*cos(phi)*sin(theta);
sin(psi)*cos(theta) cos(psi)*cos(phi) + sin(phi)*sin(theta)*sin(psi) -cos(psi)*sin(phi)+sin(theta)*sin(psi)*cos(phi);
-sin(theta) cos(theta)*sin(phi) cos(theta)*cos(phi);
];

Tnb = [
1 sin(phi)*tan(theta) cos(phi)*tan(theta);
0 cos(phi) -sin(phi);
0 sin(phi)/cos(theta) cos(phi)/cos(theta);
];

J(x,y,z,phi,theta,psi) = vpa([
Rnb zeros(3);
zeros(3) Tnb;
]);

M = MRB + MA;

C(u,v,w,p,q,r) = vpa(vpa((CRB + CA),5),5);


function y = skew(x)
    y = [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
end

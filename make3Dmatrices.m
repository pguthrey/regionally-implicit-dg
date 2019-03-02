close all
clc

syms tau real
syms xii real
syms eta real
syms zeta real

addpath('predictor')
addpath('corrector')
addpath('compute')
addpath('testfunctions')
addpath('initialize')
addpath('plotting')
addpath('numericalVN')

M = 2
r = 1
data.space_dims = 3;
data.basiscombine = 'full';
data.basis = 'canonical';
data.M = M;
data.predictor_solver = 'NewtonIteration';
data.predictorbasis = 'Q';
data.correctorbasis = 'P';

data.r_param = r; 
data.verbose = true;
data.plotIC = false; 
data.plotwhilerunning = false;
data.plotfinal = false;
data.makegif_conserved = false;
data.filter = true;%
data.smartsolver = true;

data.usewaitbars = true;
data.west_symmetry = false;
data.savefile = false;
data.email = true;
data.store = false;
data.methodtype = 'implicit';
data.Neqns = 1;

close all force

[data] = initialize_method(data);
[data] = initialize_Basis(data,data.space_dims);
[data] = initialize_gauss_quadrature(data,data.M);
%[data] = initialize_innerproduct_quadrature(data);
%[data] = initialize_mesh(data);
%[data] = initialize_precompute_testfunctions(data);
%[data] = initialize_formIntegrals(data,data.Tfinal);

matrices.theta = data.theta;
matrices.thetaT = data.thetaT;

syms xii real
testeast_term = zeros(1,data.M);
testwest_term = zeros(1,data.M);    
test = sym(zeros(1,data.M));
test_diff = sym(zeros(1,data.M));

testeast = zeros(1,data.M);
testwest = zeros(1,data.M);
for i = 1:data.M
    testeast(1,i) = testfunction_Legendre( 1,i );
    testwest(1,i) = testfunction_Legendre( -1,i );    
    test(1,i) = testfunction_Legendre(xii,i );   
    test_diff(1,i) = diff(testfunction_Legendre(xii,i ),xii); 
end
testD = zeros(data.M);
for ileft = 1:data.M
for iright = 1:data.M
   testD(ileft,iright) = double(int(test_diff(1,ileft)*test(1,iright),xii,-1,1));
end
end


[xleft,xright] = meshgrid(data.BO(:,4));
Xint = 2.*(xleft == xright);
Xeast = zeros(1,data.theta);
Xwest = zeros(1,data.theta);
for i = 1:data.thetaT
    Xeast(1,i) = testeast(data.BO(i,4));
    Xwest(1,i) = testwest(data.BO(i,4));
end
Xeasteast_term = Xeast'*Xwest;
Xwesteast_term = Xwest'*Xeast;
Xwestwest_term = Xwest'*Xwest;

Xprime_term = zeros(data.thetaT);
for ileft = 1:data.thetaT
for iright = 1:data.thetaT
   Xprime_term(ileft,iright) = testD(data.BO(ileft,4),data.BO(iright,4));
end
end

[Yleft,Yright] = meshgrid(data.BO(:,3));
Yint = 2.*(Yleft == Yright);
Yeast_term = zeros(1,data.thetaT);
Ywest_term = zeros(1,data.thetaT);
Yeast = zeros(1,data.theta);
Ywest = zeros(1,data.theta);
for i = 1:data.thetaT
    Yeast(1,i) = testeast(data.BO(i,3));
    Ywest(1,i) = testwest(data.BO(i,3));
end
Yeasteast_term = Yeast'*Ywest;
Ywesteast_term = Ywest'*Yeast;
Ywestwest_term = Ywest'*Ywest;

Yprime_term = zeros(data.thetaT);
for ileft = 1:data.thetaT
for iright = 1:data.thetaT
   Yprime_term(ileft,iright) = testD(data.BO(ileft,3),data.BO(iright,3));
end
end

[Zleft,Zright] = meshgrid(data.BO(:,2));
Zint = 2.*(Zleft == Zright);
Zeast_term = zeros(1,data.thetaT);
Zwest_term = zeros(1,data.thetaT);
Zeast = zeros(1,data.theta);
Zwest = zeros(1,data.theta);
for i = 1:data.thetaT
    Zeast(1,i) = testeast(data.BO(i,2));
    Zwest(1,i) = testwest(data.BO(i,2));
end
Zeasteast_term = Zeast'*Zwest;
Zwesteast_term = Zwest'*Zeast;
Zwestwest_term = Zwest'*Zwest;

Zprime_term = zeros(data.thetaT);
for ileft = 1:data.thetaT
for iright = 1:data.thetaT
   Zprime_term(ileft,iright) = testD(data.BO(ileft,2),data.BO(iright,2));
end
end

[Tleft,Tright] = meshgrid(data.BO(:,1));
Tint = 2.*(Tleft == Tright);
Teast = zeros(1,data.thetaT);
Twest = zeros(1,data.thetaT);
for i = 1:data.thetaT
    Teast(1,i) = testeast(data.BO(i,1));
    Twest(1,i) = testwest(data.BO(i,1));
end
Teasteast_term = Teast'*Twest;
Twesteast_term = Twest'*Teast;

Tprime_term = zeros(data.thetaT);
for ileft = 1:data.thetaT
for iright = 1:data.thetaT
   Tprime_term(ileft,iright) = testD(data.BO(ileft,1),data.BO(iright,1));
end
end

Other = Zint.*Yint.*Tint;
X_east = Xeasteast_term.*Other;
X_othr = Xwesteast_term.*Other;
X_west = Xwestwest_term.*Other;
X_der = Xprime_term.*Other;

Other = Zint.*Xint.*Tint;
Y_east = Yeasteast_term.*Other;
Y_othr = Ywesteast_term.*Other;
Y_west = Ywestwest_term.*Other;
Y_der = Yprime_term.*Other;

Other = Yint.*Xint.*Tint;
Z_east = Zeasteast_term.*Other;
Z_othr = Zwesteast_term.*Other;
Z_west = Zwestwest_term.*Other;
Z_der = Zprime_term.*Other;

Other = Yint.*Xint.*Zint;
T_east = Teasteast_term.*Other;
T_othr = Twesteast_term.*Other;
T_der = Tprime_term.*Other;

T_mat = T_east - T_der;
X_mat = X_east - X_der;
Y_mat = Y_east - Y_der;
Z_mat = Z_east - Z_der;

nux = .75;
nuy = .75;
nuz = .75;

Lxyz = T_mat + nux*(X_mat-X_west) + nuz*(Z_mat-Z_west) + nuy*(Y_mat-Y_west);
Lxy = T_mat + nux*(X_mat-X_west) + nuz*Z_mat + nuy*(Y_mat-Y_west);
Lyz = T_mat + nuz*(Z_mat-Z_west) + nux*X_mat + nuy*(Y_mat-Y_west);
Lxz = T_mat + nuz*(Z_mat-Z_west) + nuy*Y_mat + nux*(X_mat-X_west);
Lx = T_mat + nuz*Z_mat + nuy*Y_mat + nux*(X_mat-X_west);
Ly = T_mat + nux*X_mat + nuz*Z_mat + nuy*(Y_mat-Y_west);
Lz = T_mat + nux*X_mat + nuy*Y_mat + nuz*(Z_mat-Z_west);
L = T_mat + nux*X_mat + nuy*Y_mat + nuz*Z_mat; 

%{
Wim1jm1km1 = Lxyz\Tothr;
Wim1km1 = Lxz\(Tothr + nuy*Yothr*Wim1jm1km1);
Wjm1km1 = Lyz\(Tothr + nux*Xothr*Wim1jm1km1);
Wkm1 = Lz\(Tothr + nux*Xothr*Wim1km1 + nuy*Yothr*Wjm1km1);
Wim1jm1 = Lxy\(Tothr + nuz*Zothr*Wim1jm1km1);
Wim1 = Lx\(Tothr + nuy*Yothr*Wim1jm1 + nuz*Zothr*Wim1km1);
Wjm1 = Ly\(Tothr + nux*Xothr*Wim1jm1 + nuz*Zothr*Wjm1km1);
W = L\(Tothr + nux*Xothr*Wim1 + nuy*Yothr*Wjm1 + nuz*Zothr*Wkm1);
%}

P = zeros(data.thetaT,data.thetaT,2,2,2);

M1 = L\(X_mat*Lx\Y_mat + Y_mat*Ly\X_mat);
M2 = L\(Y_mat*Ly\Z_mat + Z_mat*Lz\Y_mat);
M3 = L\(X_mat*Lx\Z_mat + Z_mat*Lz\X_mat);

P(:,:,1,1,1) = L\T_mat;
P(:,:,2,1,1) = nux*L\Xothr*(Lx\Tothr);
P(:,:,1,2,1) = nuy*L\Yothr*(Ly\Tothr);
P(:,:,1,1,2) = nuz*L\Zothr*(Lz\Tothr);
P(:,:,2,2,1) = nux*nuy*M1*(Lxy\Tothr);
P(:,:,1,2,2) = nuz*nuy*M2*(Lyz\Tothr);
P(:,:,2,1,2) = nux*nuz*M3*(Lxz\Tothr);
P(:,:,2,2,2) = nux*nuy*nuz*(M1*(Lxy\Zothr)+M2*(Ly\Xothr)+M3*(Lxz\Yothr))*(Lxyz\Tothr);



[xleft,xright] = meshgrid(data.BO(:,4),data.BOs(:,3));
Xint = 2.*(xleft == xright);

Xseast = zeros(1,data.theta);
Xswest = zeros(1,data.theta);
for i = 1:data.theta
    Xseast(1,i) = testeast(data.BOs(i,3));
    Xswest(1,i) = testwest(data.BOs(i,3));
end
Xeasteast_term = Xseast'*Xwest;
Xwesteast_term = Xswest'*Xeast;

Xprime_term = zeros(data.theta,data.thetaT);
for ileft = 1:data.theta
for iright = 1:data.thetaT
   Xprime_term(ileft,iright) = testD(data.BOs(ileft,3),data.BO(iright,4));
end
end

[Yleft,Yright] = meshgrid(data.BO(:,3),data.BOs(:,2));
Yint = 2.*(Yleft == Yright);
Yseast = zeros(1,data.theta);
Yswest = zeros(1,data.theta);
for i = 1:data.theta
    Yseast(1,i) = testeast(data.BOs(i,2));
    Yswest(1,i) = testwest(data.BOs(i,2));
end
Yeasteast_term = Yseast'*Ywest;
Ywesteast_term = Yswest'*Yeast;

Yprime_term = zeros(data.theta,data.thetaT);
for ileft = 1:data.theta
for iright = 1:data.thetaT
   Yprime_term(ileft,iright) = testD(data.BOs(ileft,2),data.BO(iright,3));
end
end

[Zleft,Zright] = meshgrid(data.BO(:,2),data.BOs(:,1));
Zint = 2.*(Zleft == Zright);
Zseast = zeros(1,data.theta);
Zswest = zeros(1,data.theta);
for i = 1:data.theta
    Zseast(1,i) = testeast(data.BOs(i,1));
    Zswest(1,i) = testwest(data.BOs(i,1));
end
Zeasteast_term = Zseast'*Zwest;
Zwesteast_term = Zswest'*Zeast;

Zprime_term = zeros(data.theta,data.thetaT);
for ileft = 1:data.theta
for iright = 1:data.thetaT
   Zprime_term(ileft,iright) = testD(data.BOs(ileft,1),data.BO(iright,2));
end
end

[Tleft,Tright] = meshgrid(data.BO(:,1),1+0.*data.BOs(:,1));
Tint = 2.*(Tleft == Tright)

Other = Zint.*Yint.*Tint;
X_east = Xeasteast_term.*Other;
X_othr = Xwesteast_term.*Other;
X_der = Xprime_term.*Other;

Other = Zint.*Xint.*Tint;
Y_east = Yeasteast_term.*Other;
Y_othr = Ywesteast_term.*Other;
Y_der = Yprime_term.*Other;

Other = Yint.*Xint.*Tint;
Z_east = Zeasteast_term.*Other;
Z_othr = Zwesteast_term.*Other;
Z_der = Zprime_term.*Other;

X_mat = X_east - X_der;
Y_mat = Y_east - Y_der;
Z_mat = Z_east - Z_der;

C = zeros(data.theta,data.thetaT,2,2,2);
C(:,:,1,1,1) = nux*X_mat+nuy*Y_mat+nuz*Z_mat;
C(:,:,2,1,1) = -nux*X_othr;
C(:,:,1,2,1) = -nuy*Y_othr;
C(:,:,1,1,2) = -nuz*Z_othr;




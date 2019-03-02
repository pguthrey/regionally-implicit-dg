function [update,upindexlist,upmeshlocslist,LHS,predictor_update,corrector_update] = initialize_constantcoefficient_update_3D(data)
% written by Pierson Guthrey

syms tau real
syms xii real
syms eta real
syms zeta real

M = data.M;
r = data.r_param;

nux = data.nuv1*data.appdata.nuf;
nuy = data.nuv2*data.appdata.nug;
nuz = data.nuv3*data.appdata.nuh;

syms xii real   
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
testD = double(int(test_diff'*test,xii,-1,1));


[Xleft,Xright] = meshgrid(data.BO(:,2));
Xint = 2.*(Xleft == Xright);
[Yleft,Yright] = meshgrid(data.BO(:,3));
Yint = 2.*(Yleft == Yright);
[Zleft,Zright] = meshgrid(data.BO(:,4));
Zint = 2.*(Zleft == Zright);
[Tleft,Tright] = meshgrid(data.BO(:,1));
Tint = 2.*(Tleft == Tright);

clear Xleft Xright
clear Yleft Yright
clear Zleft Zright
clear Tleft Tright

Xeast = zeros(1,data.theta);
Xwest = zeros(1,data.theta);
for i = 1:data.thetaT
    Xeast(1,i) = testeast(data.BO(i,2));
    Xwest(1,i) = testwest(data.BO(i,2));
end
Xeasteast_term = Xeast'*Xeast;
Xwesteast_term = Xwest'*Xeast;
Xwestwest_term = Xwest'*Xwest;

Xprime_term = zeros(data.thetaT);
for ileft = 1:data.thetaT
for iright = 1:data.thetaT
   Xprime_term(ileft,iright) = testD(data.BO(ileft,2),data.BO(iright,2));
end
end

Other = Zint.*Yint.*Tint;
X_east = Xeasteast_term.*Other;
X_othr = Xwesteast_term.*Other;
X_west = Xwestwest_term.*Other;
X_der = Xprime_term.*Other;
X_mat = X_east - X_der;

clear Xwest 
clear Xeasteast_term Xwesteast_term Xwestwest_term Xprime_term 
clear X_east X_der

Yeast = zeros(1,data.theta);
Ywest = zeros(1,data.theta);
for i = 1:data.thetaT
    Yeast(1,i) = testeast(data.BO(i,3));
    Ywest(1,i) = testwest(data.BO(i,3));
end
Yeasteast_term = Yeast'*Yeast;
Ywesteast_term = Ywest'*Yeast;
Ywestwest_term = Ywest'*Ywest;

Yprime_term = zeros(data.thetaT);
for ileft = 1:data.thetaT
for iright = 1:data.thetaT
   Yprime_term(ileft,iright) = testD(data.BO(ileft,3),data.BO(iright,3));
end
end

Other = Zint.*Xint.*Tint;
Y_east = Yeasteast_term.*Other;
Y_othr = Ywesteast_term.*Other;
Y_west = Ywestwest_term.*Other;
Y_der = Yprime_term.*Other;
Y_mat = Y_east - Y_der;
clear Ywest 
clear Yeasteast_term Ywesteast_term Ywestwest_term Yprime_term
clear Y_east Y_der

Zeast = zeros(1,data.theta);
Zwest = zeros(1,data.theta);
for i = 1:data.thetaT
    Zeast(1,i) = testeast(data.BO(i,4));
    Zwest(1,i) = testwest(data.BO(i,4));
end
Zeasteast_term = Zeast'*Zeast;
Zwesteast_term = Zwest'*Zeast;
Zwestwest_term = Zwest'*Zwest;

Zprime_term = zeros(data.thetaT);
for ileft = 1:data.thetaT
for iright = 1:data.thetaT
   Zprime_term(ileft,iright) = testD(data.BO(ileft,4),data.BO(iright,4));
end
end

Other = Yint.*Xint.*Tint;
Z_east = Zeasteast_term.*Other;
Z_othr = Zwesteast_term.*Other;
Z_west = Zwestwest_term.*Other;
Z_der = Zprime_term.*Other;
Z_mat = Z_east - Z_der;

clear Zwest 
clear Zeasteast_term Zwesteast_term Zwestwest_term Zprime_term
clear Z_east Z_der

Teast = zeros(1,data.thetaT);
Twest = zeros(1,data.thetaT);
for i = 1:data.thetaT
    Teast(1,i) = testeast(data.BO(i,1));
    Twest(1,i) = testwest(data.BO(i,1));
end
Teasteast_term = Teast'*Teast;
Twesteast_term = Twest'*Teast;

Tprime_term = zeros(data.thetaT);
for ileft = 1:data.thetaT
for iright = 1:data.thetaT
   Tprime_term(ileft,iright) = testD(data.BO(ileft,1),data.BO(iright,1));
end
end

Other = Yint.*Xint.*Zint;
T_east = Teasteast_term.*Other;
T_othr = Twesteast_term.*Other;
T_der = Tprime_term.*Other;
T_mat = T_east - T_der;

clear Teast Twest 
clear Teasteast_term Twesteast_term Tprime_term
clear T_east T_der
clear Xint Yint Zint Tint
disp('matrices crafted')

invLxyz = inv(T_mat + nux*(X_mat-X_west) + nuz*(Z_mat-Z_west) + nuy*(Y_mat-Y_west));
switch data.r_param
    case 1
        %{
        Lxyz = (T_mat + nux*(X_mat-X_west) + nuz*(Z_mat-Z_west) + nuy*(Y_mat-Y_west));
        Lxy = (T_mat + nux*(X_mat-X_west) + nuz*Z_mat + nuy*(Y_mat-Y_west));
        Lyz = (T_mat + nuz*(Z_mat-Z_west) + nux*X_mat + nuy*(Y_mat-Y_west));
        Lxz = (T_mat + nuz*(Z_mat-Z_west) + nuy*Y_mat + nux*(X_mat-X_west));
        Lx = (T_mat + nuz*Z_mat + nuy*Y_mat + nux*(X_mat-X_west));
        Ly = (T_mat + nux*X_mat + nuz*Z_mat + nuy*(Y_mat-Y_west));
        Lz = (T_mat + nux*X_mat + nuy*Y_mat + nuz*(Z_mat-Z_west));
        L = (T_mat + nux*X_mat + nuy*Y_mat + nuz*Z_mat); 
        
        X = nux*X_othr;
        Y = nuy*Y_othr;
        Z = nuz*Z_othr;
        L0 = Lxyz.*0;
        
        %}
        LHS = 1;
        %{
        [[Lxyz L0 L0 L0 L0 L0 L0 L0];
                [-X Lyz L0 L0 L0 L0 L0 L0];
                [-Y L0 Lxz L0 L0 L0 L0 L0];
                [L0 -Y -X Lz L0 L0 L0 L0];
                [-Z L0 L0 L0 Lxy L0 L0 L0];
                [L0 -Z L0 L0 -X Ly L0 L0];
                [L0 L0 -Z L0 -Y L0 Lx L0];
                [L0 L0 L0 -Z L0 -Y -X L];]
                ;
        %}
                     
        invLxy = inv(T_mat + nux*(X_mat-X_west) + nuz*Z_mat + nuy*(Y_mat-Y_west));
        invLyz = inv(T_mat + nuz*(Z_mat-Z_west) + nux*X_mat + nuy*(Y_mat-Y_west));
        invLxz = inv(T_mat + nuz*(Z_mat-Z_west) + nuy*Y_mat + nux*(X_mat-X_west));
        invLx = inv(T_mat + nuz*Z_mat + nuy*Y_mat + nux*(X_mat-X_west));
        invLy = inv(T_mat + nux*X_mat + nuz*Z_mat + nuy*(Y_mat-Y_west));
        invLz = inv(T_mat + nux*X_mat + nuy*Y_mat + nuz*(Z_mat-Z_west));
        invL = inv(T_mat + nux*X_mat + nuy*Y_mat + nuz*Z_mat); 

        clear T_mat
        clear X_mat Y_mat Z_mat
        clear X_west Y_west Z_west

        
        predictor_update = zeros(data.thetaT,data.thetaT,2,2,2);


        
        disp('matrices inverted')
        
        predictor_update(:,:,1,1,2) = nuz*(invL*(Z_othr*(invLz)));
        predictor_update(:,:,1,2,1) = nuy*(invL*(Y_othr*(invLy)));
        predictor_update(:,:,2,1,1) = nux*(invL*(X_othr*(invLx)));    
        M1 = invL*(X_othr*(invLx*Y_othr) + Y_othr*(invLy*X_othr));
        M2 = invL*(Y_othr*(invLy*Z_othr) + Z_othr*(invLz*Y_othr));
        clear invLy 
        M3 = invL*(X_othr*(invLx*Z_othr) + Z_othr*(invLz*X_othr));
        clear invLx
        clear invLz
        predictor_update(:,:,2,2,2) = nux*nuy*nuz*(M1*(invLxy)*Z_othr+M2*(invLyz)*X_othr+M3*(invLxz)*Y_othr)*(invLxyz);
        clear Lxyz
        predictor_update(:,:,2,2,1) = nux*nuy*M1*(invLxy);
        clear M1 invLxy
        predictor_update(:,:,1,2,2) = nuz*nuy*M2*(invLyz);
        clear M2 invLyz
        predictor_update(:,:,2,1,2) = nux*nuz*M3*(invLxz);
        clear M3 invLxz

        clear invLx X_othr

        predictor_update(:,:,1,1,1) = invL;
                
        %clear M1 M2 M3
        clear L 
        clear Lxy 
        %clear Lxyz Lyz Lxy Lz

        
    case 0
        predictor_update = zeros(data.thetaT,data.thetaT,1,1,1);
        predictor_update(:,:,1,1,1) = invLxyz;
        LHS = inv(invLxyz);
        
end

disp('making corrector')


[Tleft,Tright] = meshgrid(data.BO(:,1),1+0.*data.BOs(:,1));
Tint = 2.*(Tleft == Tright);
[Xleft,Xright] = meshgrid(data.BO(:,2),data.BOs(:,1));
Xint = 2.*(Xleft == Xright);
[Yleft,Yright] = meshgrid(data.BO(:,3),data.BOs(:,2));
Yint = 2.*(Yleft == Yright);
[Zleft,Zright] = meshgrid(data.BO(:,4),data.BOs(:,3));
Zint = 2.*(Zleft == Zright);

Xseast = zeros(1,data.theta);
Xswest = zeros(1,data.theta);
for i = 1:data.theta
    Xseast(1,i) = testeast(data.BOs(i,1));
    Xswest(1,i) = testwest(data.BOs(i,1));
end
Xeasteast_term = Xseast'*Xeast;
Xwesteast_term = Xswest'*Xeast;

Xprime_term = zeros(data.theta,data.thetaT);
for ileft = 1:data.theta
for iright = 1:data.thetaT
   Xprime_term(ileft,iright) = testD(data.BOs(ileft,1),data.BO(iright,2));
end
end

Other = Zint.*Yint.*Tint;
X_east = Xeasteast_term.*Other;
X_othr = Xwesteast_term.*Other;
X_der = Xprime_term.*Other;
X_mat = X_east - X_der;

clear Xswest Xseast Xeast
clear Xeasteast_term Xwesteast_term Xwestwest_term Xprime_term 
clear X_east X_der


Yseast = zeros(1,data.theta);
Yswest = zeros(1,data.theta);
for i = 1:data.theta
    Yseast(1,i) = testeast(data.BOs(i,2));
    Yswest(1,i) = testwest(data.BOs(i,2));
end
Yeasteast_term = Yseast'*Yeast;
Ywesteast_term = Yswest'*Yeast;

Yprime_term = zeros(data.theta,data.thetaT);
for ileft = 1:data.theta
for iright = 1:data.thetaT
   Yprime_term(ileft,iright) = testD(data.BOs(ileft,2),data.BO(iright,3));
end
end

Other = Zint.*Xint.*Tint;
Y_east = Yeasteast_term.*Other;
Y_othr = Ywesteast_term.*Other;
Y_der = Yprime_term.*Other;
Y_mat = Y_east - Y_der;

clear Yswest Yseast Yeast
clear Yeasteast_term Ywesteast_term Ywestwest_term Yprime_term 
clear Y_east Y_der

Zseast = zeros(1,data.theta);
Zswest = zeros(1,data.theta);
for i = 1:data.theta
    Zseast(1,i) = testeast(data.BOs(i,3));
    Zswest(1,i) = testwest(data.BOs(i,3));
end
Zeasteast_term = Zseast'*Zeast;
Zwesteast_term = Zswest'*Zeast;

Zprime_term = zeros(data.theta,data.thetaT);
for ileft = 1:data.theta
for iright = 1:data.thetaT
   Zprime_term(ileft,iright) = testD(data.BOs(ileft,3),data.BO(iright,4));
end
end

Other = Yint.*Xint.*Tint;
Z_east = Zeasteast_term.*Other;
Z_othr = Zwesteast_term.*Other;
Z_der = Zprime_term.*Other;
Z_mat = Z_east - Z_der;

clear Zswest Zseast Zeast
clear Zeasteast_term Zwesteast_term Zwestwest_term Zprime_term 
clear Z_east Z_der

corrector_update = zeros(data.theta,data.thetaT,2,2,2);
corrector_update(:,:,1,1,1) = nux*X_mat+nuy*Y_mat+nuz*Z_mat;
corrector_update(:,:,2,1,1) = -nux*X_othr;
corrector_update(:,:,1,2,1) = -nuy*Y_othr;
corrector_update(:,:,1,1,2) = -nuz*Z_othr;

corrector_update = corrector_update/2^(data.space_dims);

[Tprojl,Tprojr] = meshgrid(1+0.*data.BOs(:,1),data.BO(:,1));
[Xprojl,Xprojr] = meshgrid(data.BOs(:,1),data.BO(:,2));
[Yprojl,Yprojr] = meshgrid(data.BOs(:,2),data.BO(:,3));
[Zprojl,Zprojr] = meshgrid(data.BOs(:,3),data.BO(:,4));
Proj = (Tprojl == Tprojr);
Proj = Proj.*(Xprojl == Xprojr);
Proj = Proj.*(Yprojl == Yprojr);
Proj = Proj.*(Zprojl == Zprojr);
clear Tprojl Tprojr
clear Xprojl Xprojr
clear Yprojl Yprojr
clear Zprojl Zprojr


disp('corrector matrices made')


upindexlist = [1 1 1];
upmeshlocslist = [0 0 0];
update = zeros(data.theta,data.theta,r+2,r+2,r+2);

crrlist = [ [0 0 0];
            [-1 0 0];
            [0 -1 0];
            [0 0 -1]];

for ic = 1:length(crrlist(:,1))
    icx = crrlist(ic,1);
    icy = crrlist(ic,2);
    icz = crrlist(ic,3);
    C = corrector_update(:,:,1+abs(icx),1+abs(icy),1+abs(icz));
    for irz = -r:0
    for iry = -r:0
    for irx = -r:0
        P = predictor_update(:,:,1+abs(irx),1+abs(iry),1+abs(irz));        
        iux = 1 + abs(irx) + abs(icx);
        iuy = 1 + abs(iry) + abs(icy);
        iuz = 1 + abs(irz) + abs(icz); 
        mat = C*P*T_othr*Proj;
        if norm(mat) > 0
            upindex = [iux iuy iuz];
            check1 = upindexlist(:,1)==upindex(1);
            check2 = upindexlist(:,2)==upindex(2);
            check3 = upindexlist(:,3)==upindex(3);
            check = sum(check1.*check2.*check3);
            if ~check
                upindexlist = [upindexlist ; upindex];
                upmeshlocs = [icx+irx icy+iry icz+irz];
                upmeshlocslist = [upmeshlocslist ; upmeshlocs ];
            end
            update(:,:,iux,iuy,iuz) = update(:,:,iux,iuy,iuz) + mat;
        end
    end
    end
    end
end
disp('update made')




end


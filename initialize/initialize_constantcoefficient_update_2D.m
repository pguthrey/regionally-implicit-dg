function [update,upindexlist,upmeshlocslist,LHS,predictor_update,corrector_update,mats] = initialize_constantcoefficient_update_2D(data)
% written by Pierson Guthrey


M = data.M;
r = data.r_param;

nux = data.nuv1*data.appdata.nuf;
nuy = data.nuv2*data.appdata.nug;

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

[xleft,xright] = meshgrid(data.BO(:,2));
Xint = 2.*(xleft == xright);
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

[Yleft,Yright] = meshgrid(data.BO(:,3));
Yint = 2.*(Yleft == Yright);
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

[Tleft,Tright] = meshgrid(data.BO(:,1));
Tint = 2.*(Tleft == Tright);
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

Other = Yint.*Tint;
X_east = Xeasteast_term.*Other;
X_othr = Xwesteast_term.*Other;
X_west = Xwestwest_term.*Other;
X_der = Xprime_term.*Other;

Other = Xint.*Tint;
Y_east = Yeasteast_term.*Other;
Y_othr = Ywesteast_term.*Other;
Y_west = Ywestwest_term.*Other;
Y_der = Yprime_term.*Other;

Other = Yint.*Xint;
T_east = Teasteast_term.*Other;
T_othr = Twesteast_term.*Other;
T_der = Tprime_term.*Other;

T_mat = T_east - T_der;
X_mat = X_east - X_der;
Y_mat = Y_east - Y_der;


Lxy = T_mat + nux*(X_mat-X_west) + nuy*(Y_mat-Y_west);
Lx = T_mat + nuy*Y_mat + nux*(X_mat-X_west);
Ly = T_mat + nux*X_mat + nuy*(Y_mat-Y_west);
L = T_mat + nux*X_mat + nuy*Y_mat; 


switch data.r_param
    case 1
        Lxy = T_mat + nux*(X_mat-X_west) + nuy*(Y_mat-Y_west);
        Lx = T_mat + nuy*Y_mat + nux*(X_mat-X_west);
        Ly = T_mat + nux*X_mat + nuy*(Y_mat-Y_west);
        L = T_mat + nux*X_mat + nuy*Y_mat; 

        P = zeros(data.thetaT,data.thetaT,2,2,1);

        M1 = L\(X_othr*(Lx\Y_othr) + Y_othr*(Ly\X_othr));

        predictor_update(:,:,1,1,1) = inv(L);
        predictor_update(:,:,2,1,1) = nux*(L\(X_othr*inv(Lx)));
        predictor_update(:,:,1,2,1) = nuy*(L\(Y_othr*inv(Ly)));
        predictor_update(:,:,2,2,1) = nux*nuy*M1*inv(Lxy);
  
        LHS = [[Lxy 0.*L 0.*L 0.*L]
               [-nux*X_othr Ly 0.*L 0.*L]
               [-nuy*Y_othr 0.*L Lx 0.*L]
               [0.*L -nuy*Y_othr -nux*X_othr L]];
    case 0
        LHS = Lxy;
        P = zeros(data.thetaT,data.thetaT,1,1,1);
        predictor_update(:,:,1,1,1) = inv(Lxy);
end

Z = 0.*T_mat;
LHS_time = [[T_mat Z Z Z];[Z T_mat Z Z];[Z Z T_mat Z];[Z Z Z T_mat];] ;
LHS_x = [[X_mat-X_west Z Z Z];[-X_othr X_mat Z Z];[Z Z X_mat-X_west Z];[Z Z -X_othr X_mat];] ;
LHS_y = [[Y_mat-Y_west Z Z Z];[Z Y_mat-Y_west Z Z];[-Y_othr Z Y_mat Z];[Z -Y_othr Z Y_mat];] ;


LHS_tilde_t = LHS_time\LHS_time;
LHS_tilde_x = inv(LHS_time)*LHS_x;
LHS_tilde_y = inv(LHS_time)*LHS_y;

clc
close all
SO = 30
nu = data.nuv1
neumann = zeros(size(T_mat));
Itilde = T_mat - nu*(X_east+Y_east)
Ttilde = nu*(X_der+Y_der) 
exact = inv(Itilde-Ttilde);
T = inv(Itilde)*Ttilde
for M=0:SO
    neumann = neumann + T^M;
    error(1+M) = norm(neumann*inv(Itilde)-exact)
end
semilogy(0:SO,error)



keyboard 

SO = 8
nu =data.nuv1
exact = inv(LHS_time+nu*(LHS_x+LHS_y));
neumann = zeros(size(exact));
T = -nu*(LHS_tilde_x+LHS_tilde_y)
for M=0:SO
    neumann = neumann + T^M;
    error(1+M) = norm(neumann*inv(LHS_time)-exact)
end
semilogy(0:SO,error)
    
LHSinv = neumann*inv(LHS_time);
indices = 1:data.thetaT;
predictor_update(:,:,1,1,1) = LHSinv(3*data.thetaT+indices,3*data.thetaT+indices);
predictor_update(:,:,2,1,1) = LHSinv(3*data.thetaT+indices,2*data.thetaT+indices);
predictor_update(:,:,1,2,1) = LHSinv(3*data.thetaT+indices,data.thetaT+indices);
predictor_update(:,:,2,2,1) = LHSinv(3*data.thetaT+indices,indices);


[xleft,xright] = meshgrid(data.BO(:,2),data.BOs(:,1));
Xint = 2.*(xleft == xright);

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

[Yleft,Yright] = meshgrid(data.BO(:,3),data.BOs(:,2));
Yint = 2.*(Yleft == Yright);
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

[Tleft,Tright] = meshgrid(data.BO(:,1),1+0.*data.BOs(:,1));
Tint = 2.*(Tleft == Tright);

Other = Yint.*Tint;
X_east = Xeasteast_term.*Other;
X_othr = Xwesteast_term.*Other;
X_der = Xprime_term.*Other;

Other = Xint.*Tint;
Y_east = Yeasteast_term.*Other;
Y_othr = Ywesteast_term.*Other;
Y_der = Yprime_term.*Other;

X_mat = X_east - X_der;
Y_mat = Y_east - Y_der;

corrector_update = zeros(data.theta,data.thetaT,2,2);
corrector_update(:,:,1,1,1) = nux*X_mat+nuy*Y_mat;
corrector_update(:,:,2,1,1) = -nux*X_othr;
corrector_update(:,:,1,2,1) = -nuy*Y_othr;

corrector_update = corrector_update/2^(data.space_dims);

[Tprojl,Tprojr] = meshgrid(1+0.*data.BOs(:,1),data.BO(:,1));
[Xprojl,Xprojr] = meshgrid(data.BOs(:,1),data.BO(:,2));
[Yprojl,Yprojr] = meshgrid(data.BOs(:,2),data.BO(:,3));
Proj = (Tprojl == Tprojr);
Proj = Proj.*(Xprojl == Xprojr);
Proj = Proj.*(Yprojl == Yprojr);

upindexlist = [1 1 1];
upmeshlocslist = [0 0 0];
update = zeros(data.theta,data.theta,r+2,r+2);

crrlist = [ [0 -1 0];
            [-1 0 0];
            [0 0 0];];
mats.Umatlist2 = [];
mats.Pmatlist2 = [];
mats.Cmatlist2 = [];
mats.Tmatlist2 = [];
for ic = 1:length(crrlist(:,1))
    icx = crrlist(ic,1);
    icy = crrlist(ic,2);
    icz = crrlist(ic,3);
    C = corrector_update(:,:,1+abs(icx),1+abs(icy),1+abs(icz));
    for irz = 0
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
            check = sum(check1.*check2);            
            if ~check
                upindexlist = [upindexlist ; upindex];
                upmeshlocs = [icx+irx icy+iry icz+irz];
                upmeshlocslist = [upmeshlocslist ; upmeshlocs ];
            end
            mats.Pmatlist2 = [mats.Pmatlist2 ; P ];
            mats.Cmatlist2 = [mats.Cmatlist2 ; C ];
            mats.Umatlist2 = [mats.Umatlist2 ; mat ];
            mats.Tmatlist2 = [mats.Tmatlist2 ; T_othr*Proj ];
            update(:,:,iux,iuy,iuz) = update(:,:,iux,iuy,iuz) + mat;
        end
    end
    end
    end
end


end


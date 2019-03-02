function [update,upindexlist,upmeshlocslist,LHS,predictor_update,corrector_update] = initialize_constantcoefficient_update_1D(data)
% written by Pierson Guthrey

M = data.M;
r = data.r_param;

data.appdata.nuf = 1;
nux = data.nuv1*data.appdata.nuf;

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

Other = Tint;
X_east = Xeasteast_term.*Other;
X_othr = Xwesteast_term.*Other;
X_west = Xwestwest_term.*Other;
X_der = Xprime_term.*Other;

Other = Xint;
T_east = Teasteast_term.*Other;
T_othr = Twesteast_term.*Other;
T_der = Tprime_term.*Other;

T_mat = T_east - T_der;
X_mat = X_east - X_der;

Lx = T_mat + nux*(X_mat-X_west);

switch data.r_param
    case 1
        Lx = T_mat + nux*(X_mat-X_west);
        L = T_mat + nux*X_mat; 

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

        P = sym(zeros(data.thetaT,data.thetaT,2,1,1));

        predictor_update(:,:,1,1,1) = inv(L);
        predictor_update(:,:,2,1,1) = nux*(L\(X_othr*inv(Lx)));
  
        LHS =  [[Lx 0.*L]
               [-nux*X_othr L]];
    case 0
        LHS = Lx;
        P = zeros(data.thetaT,data.thetaT,1,1,1);
        predictor_update(:,:,1,1,1) = inv(Lx);
end

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

[Tleft,Tright] = meshgrid(data.BO(:,1),1+0.*data.BOs(:,1));
Tint = 2.*(Tleft == Tright);

Other = Tint;
X_east = Xeasteast_term.*Other;
X_othr = Xwesteast_term.*Other;
X_der = Xprime_term.*Other;

X_mat = X_east - X_der;

corrector_update = sym(zeros(data.theta,data.thetaT,2));
corrector_update(:,:,1,1,1) = nux*X_mat;
corrector_update(:,:,2,1,1) = -nux*X_othr;

corrector_update = corrector_update/2^(data.space_dims);

[Tprojl,Tprojr] = meshgrid(1+0.*data.BOs(:,1),data.BO(:,1));
[Xprojl,Xprojr] = meshgrid(data.BOs(:,1),data.BO(:,2));
Proj = (Tprojl == Tprojr);
Proj = Proj.*(Xprojl == Xprojr);

upindexlist = [1 1 1];
upmeshlocslist = [0 0 0];
update = sym(zeros(data.theta,data.theta,r+2));

crrlist = [ [0 0 0];
            [-1 0 0];]

for ic = 1:length(crrlist(:,1))
    icx = crrlist(ic,1);
    icy = crrlist(ic,2);
    icz = crrlist(ic,3);
    C = corrector_update(:,:,1+abs(icx),1+abs(icy),1+abs(icz));
    for irz = 0
    for iry = 0
    for irx = -r:0
        P = predictor_update(:,:,1+abs(irx),1+abs(iry),1+abs(irz));        
        iux = 1 + abs(irx) + abs(icx);
        iuy = 1 + abs(iry) + abs(icy);
        iuz = 1 + abs(irz) + abs(icz); 
        mat = C*P*T_othr*Proj;
        if 1%norm(mat) > 0
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

keyboard 

end


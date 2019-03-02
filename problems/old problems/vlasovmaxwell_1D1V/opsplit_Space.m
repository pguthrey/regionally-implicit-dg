function [fcoeffs_new,E1coeffs_new] = opsplit_Space(deltat,fcoeffs,E1coeffs,data)

h1 = waitbar(0,'Advecting in Space','OuterPosition', [800 100 300 75]);

gamma = @(v1,v2) 1;

normphi = 4;
fcoeffs_new = zeros(data.theta1D1V,data.Nx,data.Nv1);

for iv1 = 1:data.Nv1
    for v1quad = 1:data.P
        eta1 = data.locs(v1quad);
        v1 = data.v1centers(iv1)+eta1*data.deltav1/2;
        
        LeftMat = zeros(data.theta1D1V,data.theta1D1V);
        RiteMat = zeros(data.theta1D1V,data.theta1D1V);
        
        deltas = v1*deltat/gamma(v1);
        jay = floor(deltas/data.deltax);
        nu = deltas/data.deltax - jay;
        
        for xquad = 1:data.P
            xii = data.locs(xquad);
            xiitilde1 = (1-nu)*xii+nu;
            xiitilde2 = (1-nu)*xii-nu;
            xiitilde3 = nu*xii+nu-1;
            xiitilde4 = nu*xii-nu+1;
            dxii1 = (1-nu)*data.wgts2D(xquad,v1quad);
            dxii2 = nu*data.wgts2D(xquad,v1quad);
        for ir = 1:data.theta1D1V
        for ic = 1:data.theta1D1V
            phil = testfunction_zeta1D1V(xiitilde1,eta1,ir,data);
            phik = testfunction_zeta1D1V(xiitilde2,eta1,ic,data);
            RiteMat(ir,ic) = RiteMat(ir,ic) + phil*phik*dxii1;
            phil = testfunction_zeta1D1V(xiitilde3,eta1,ir,data);
            phik = testfunction_zeta1D1V(xiitilde4,eta1,ic,data);
            LeftMat(ir,ic) = LeftMat(ir,ic) + phil*phik*dxii2;
        end
        end
        end        
        Fimjm1 = circshift(fcoeffs,jay+1,2);
        Fimj = circshift(fcoeffs,jay,2);        
        for ix = 1:data.Nx
            Fnew = (LeftMat*Fimjm1(:,ix,iv1)+RiteMat*Fimj(:,ix,iv1))/normphi;
            fcoeffs_new(:,ix,iv1) = fcoeffs_new(:,ix,iv1) + Fnew;
        end
    end
end

%{
J1old = @(x,v1) v1*fnaut(x,v1);
J1coeffs_old = project_1D1V_to_1D(J1old,data);

%J1new = @(x,v1) v1*fnaut(x-deltat*v1/gamma(v1),v1);
%J1coeffs_new = project_1D1V_to_1D(J1new,data);

J1coeffs_new =zeros(data.theta1D,data.Nx);


normphi = 2;
for iv1 = 1:data.Nv1
    for v1quad = 1:data.P
        eta1 = data.locs(v1quad);
        v1 = data.v1centers(iv1)+eta1*data.deltav1/2;
        
        LeftMat = zeros(data.theta1D,data.theta1D1V);
        RiteMat = zeros(data.theta1D,data.theta1D1V);
        
        deltas = v1*deltat/gamma(v1);
        jay = floor(deltas/data.deltax);
        nu = deltas/data.deltax - jay;
        
        for xquad = 1:data.P
            xii = data.locs(xquad);
            xiitilde1 = nu*xii+1-nu;
            xiitilde2 = nu*xii+nu-1;
            xiitilde3 = (1-nu)*xii-nu;
            xiitilde4 = (1-nu)*xii+nu;
            dxii1 = nu*data.wgts2D(xquad,v1quad);
            dxii2 = (1-nu)*data.wgts2D(xquad,v1quad);
        for ir = 1:data.theta1D;
        for ic = 1:data.theta1D1V;
            phil = testfunction_zeta1D(xiitilde1,ir,data);
            phik = testfunction_zeta1D1V(xiitilde2,eta1,ic,data);
            LeftMat(ir,ic) = LeftMat(ir,ic) + phil*phik*dxii1;
            phil = testfunction_zeta1D(xiitilde3,ir,data);
            phik = testfunction_zeta1D1V(xiitilde4,eta1,ic,data);
            RiteMat(ir,ic) = RiteMat(ir,ic) + phil*phik*dxii2;
        end
        end
        end
        
        Fimjm1 = circshift(fcoeffs,jay+1,2);
        Fimj = circshift(fcoeffs,jay,2);
        
        for ix = 1:data.Nx
            Fnew = v1*(LeftMat*Fimjm1(:,ix,iv1)+RiteMat*Fimj(:,ix,iv1))/normphi;
            J1coeffs_new(:,ix) = J1coeffs_new(:,ix) + Fnew;
        end
    end
end
%}

%data = VM1D1V_makegif(fcoeffs,0.*E1coeffs,0,data,0);
%data = VM1D1V_makegif(fcoeffs_new,0.*E1coeffs,0,data,0);

%keyboard 

J1coeffs_new = zeros(data.theta1D,data.Nx);
normphi = 2;
for ix = 1:data.Nx
    temp = zeros(data.theta1D,1);
    for iv1 = 1:data.Nv1
        for v1quad = 1:data.P
        for xquad = 1:data.P
            eta1 = data.locs(v1quad);
            %x = data.xcenters(ix)+xii*data.deltax/2;
            v1 = data.v1centers(data.Nv1+1-iv1)+eta1*data.deltav1/2;
            zetatrans = data.zeta1D(1,:,xquad)';
            zeta1D1V = data.zeta1D1V(1,:,xquad,v1quad);
            wgt = data.wgts2D(xquad,v1quad);
            coeffs = fcoeffs_new(:,ix,iv1);
            temp = temp + v1*zetatrans*zeta1D1V*coeffs*wgt/normphi;
        end
        end
    end
    J1coeffs_new(:,ix) = temp;
    %waitbar(ix/Nx,h1)
end

J1coeffs_old = zeros(data.theta1D,data.Nx);
normphi = 2;
for ix = 1:data.Nx
    temp = zeros(data.theta1D,1);
    for iv1 = 1:data.Nv1
        for v1quad = 1:data.P
        for xquad = 1:data.P
            eta1 = data.locs(v1quad);
            %x = data.xcenters(ix)+xii*data.deltax/2;
            v1 = data.v1centers(data.Nv1+1-iv1)+eta1*data.deltav1/2;
            zetatrans = data.zeta1D(1,:,xquad)';
            zeta1D1V = data.zeta1D1V(1,:,xquad,v1quad);
            wgt = data.wgts2D(xquad,v1quad);
            coeffs = fcoeffs(:,ix,iv1);
            temp = temp + v1*zetatrans*zeta1D1V*coeffs*wgt/normphi;
        end
        end
    end
    J1coeffs_old(:,ix) = temp;
    %waitbar(ix/Nx,h1)
end

%Euler Imp
E1coeffs_new = E1coeffs-(J1coeffs_new)*deltat;
%Trapezoid Rule
%E1coeffs_new = E1coeffs-(J1coeffs_old+J1coeffs_new)*deltat/2;


%{
close all
plot(J1coeffs_new(1,:))
figure
plot(J1coeffs_new(2,:))
figure
plot(J1coeffs_old(1,:))
figure
plot(J1coeffs_old(2,:))
figure
plot(E1coeffs(1,:))
figure
plot(E1coeffs(2,:)) 

figure
plot(E1coeffs_new(1,:))
figure
plot(E1coeffs_new(2,:))

keyboard
%}

%keyboard 

delete(h1)

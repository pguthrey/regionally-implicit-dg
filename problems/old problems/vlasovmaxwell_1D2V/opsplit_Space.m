function [fcoeffs_new,E1coeffs_new,A1coeffs_new,A2coeffs_new] = opsplit_Space(deltat,fcoeffs,E1coeffs,A1coeffs,A2coeffs,data)

h1 = waitbar(0,'Advecting in Space','OuterPosition', [800 100 300 75]);

fnaut = @(x,v1,v2) DGeval_1D2V(fcoeffs,x,v1,v2,data); 
A1naut = @(x) DGeval_1D(A1coeffs,x,data);
A2naut = @(x) DGeval_1D(A2coeffs,x,data);

qe = 1; %
Rm = 1;

gamma = @(v1,v2) 1;

normphi = 8;
fcoeffs_new = zeros(data.theta1D2V,data.Nx,data.Nv1,data.Nv2);

for iv1 = 1:data.Nv1
for iv2 = 1:data.Nv2
    for v1quad = 1:data.P
    for v2quad = 1:data.P
        eta1 = data.locs(v1quad);
        eta2 = data.locs(v2quad);
        v1 = data.v1centers(iv1)+eta1*data.deltav1/2;
        v2 = data.v2centers(iv2)+eta2*data.deltav2/2;
        
        LeftMat = zeros(data.theta1D2V,data.theta1D2V);
        RiteMat = zeros(data.theta1D2V,data.theta1D2V);
        
        deltas = v1*deltat/gamma(v1,v2);
        jay = floor(deltas/data.deltax);
        nu = deltas/data.deltax - jay;
        
        for xquad = 1:data.P
            xii = data.locs(xquad);
            xiitilde1 = nu*xii+1-nu;
            xiitilde2 = nu*xii+nu-1;
            xiitilde3 = (1-nu)*xii-nu;
            xiitilde4 = (1-nu)*xii+nu;
            dxii1 = nu*data.wgts3D(xquad,v1quad,v2quad);
            dxii2 = (1-nu)*data.wgts3D(xquad,v1quad,v2quad);
        for ir = 1:data.theta1D2V;
        for ic = 1:data.theta1D2V;
            phil = testfunction_zeta1D2V(xiitilde1,eta1,eta2,ir,data);
            phik = testfunction_zeta1D2V(xiitilde2,eta1,eta2,ic,data);
            LeftMat(ir,ic) = LeftMat(ir,ic) + phil*phik*dxii1;
            phil = testfunction_zeta1D2V(xiitilde3,eta1,eta2,ir,data);
            phik = testfunction_zeta1D2V(xiitilde4,eta1,eta2,ic,data);
            RiteMat(ir,ic) = RiteMat(ir,ic) + phil*phik*dxii2;
        end
        end
        end
        
        Fimjm1 = circshift(fcoeffs,jay+1,2);
        Fimj = circshift(fcoeffs,jay,2);
        
        for ix = 1:data.Nx
            Fnew = (LeftMat*Fimjm1(:,ix,iv1,iv2)+RiteMat*Fimj(:,ix,iv1,iv2))/normphi;
            fcoeffs_new(:,ix,iv1,iv2) = fcoeffs_new(:,ix,iv1,iv2) + Fnew;
        end
                
    end
    end
end
end

%fadvec = @(x,v1,v2) DGeval_1D2V(fcoeffs,x-v1.*deltat/gamma(v1,v2),v1,v2,data); 
%fcoeffs_new = project_1D2V_to_1D2V(fadvec,data);

J1old = @(x,v1,v2) v1*fnaut(x,v1,v2);
J1new = @(x,v1,v2) v1*fnaut(x-deltat*v1/gamma(v1,v2),v1,v2);

J1coeffs_old = project_1D2V_to_1D(J1old,data);
J1coeffs_new = project_1D2V_to_1D(J1new,data);

J2old = @(x,v1,v2) v2*fnaut(x,v1,v2);
J2new = @(x,v1,v2) v2*fnaut(x-deltat*v1/gamma(v1,v2),v1,v2);

J2coeffs_old = project_1D2V_to_1D(J2old,data);
J2coeffs_new = project_1D2V_to_1D(J2new,data);

A1advec = project_1D_to_1D(@(x) A1naut(x+deltat),data);
A2advec = project_1D_to_1D(@(x) A2naut(x-deltat),data);

J2newadvec = @(x) DGeval_1D(J2coeffs_new,x,data);

J2newadvecleft = project_1D_to_1D(@(x) J2newadvec(x+deltat),data);
J2newadvecrite = project_1D_to_1D(@(x) J2newadvec(x-deltat),data);

integrateJ2left = (J2newadvecleft+J2coeffs_old)*deltat/2;
integrateJ2rite = (J2newadvecrite+J2coeffs_old)*deltat/2;

E1coeffs_new = E1coeffs-(J1coeffs_old+J1coeffs_new)*deltat/2;
A1coeffs_new = A1advec + integrateJ2left;
A2coeffs_new = A2advec - integrateJ2rite;

delete(h1)

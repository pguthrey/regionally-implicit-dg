function [Jac_cell,residual_cell,maxspeedF,maxspeedG] = res_and_Jac_brute(DGregion,DGregion_past,nuv1,nuv2,data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

vectpsi = data.vectpsi;
vectpsi_dxii_Trans = data.vectpsi_dxii_Trans;
vectpsi_deta_Trans = data.vectpsi_deta_Trans;


tauflux = data.tauflux;
Ipast = data.I_past;
thetaT = data.thetaT;
Pdp1 = data.Pdp1;
space_dims = data.space_dims;
Dp1list = data.Dp1list;
Dp1quadwgts = data.Dp1quadwgts;

residual_cell = zeros(thetaT,3,3);
Jac_cell = zeros(thetaT,thetaT,3,3);

v3quad = 1;
maxspeedF = 0;
maxspeedG = 0;

for rx = 1:3
    for ry = 1:3
        qstar = DGregion(:,rx,ry);
        residual_cell(:,rx,ry) = tauflux*qstar-Ipast*DGregion_past(:,rx,ry);
        Jac_cell(:,:,rx,ry) = tauflux; 
    end
end
        
for k = 1:Pdp1
        tquad = Dp1list(k,1); 
        v1quad = Dp1list(k,2); 
        if space_dims >= 2
            v2quad = Dp1list(k,3);
            if space_dims >= 3
                v3quad = Dp1list(k,4);
            end
        end
        weight = Dp1quadwgts(k);

        psidxii = vectpsi_dxii_Trans(:,:,tquad,v1quad,v2quad,v3quad);
        psideta = vectpsi_deta_Trans(:,:,tquad,v1quad,v2quad,v3quad);
        psi = vectpsi(:,:,tquad,v1quad,v2quad,v3quad);
        
        mult = nuv1*weight*psidxii + nuv2*weight*psideta;
        %psidxiipsi = psidxii*psi;
        %psidetapsi = psideta*psi;
        
        for rx = 1:3
            for ry = 1:3
                %Derivative terms
                qstar = DGregion(:,rx,ry);
                wstar = psi*qstar;
                Fstar = wstar^2/2;
                term = mult*Fstar;
                                
                residual_cell(:,rx,ry) = residual_cell(:,rx,ry) - term; 
                maxspeedF = max(maxspeedF,abs(wstar));
                maxspeedG = max(maxspeedG,abs(wstar));

                term = mult*wstar*psi;
                Jac_cell(:,:,rx,ry) = Jac_cell(:,:,rx,ry) - term;
            end
        end
end 



end


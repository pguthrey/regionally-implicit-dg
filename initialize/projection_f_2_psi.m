function [ DGcoeffs ] = projection_f_2_psi(f,tnow,deltat,data)
% Projects f onto the spacetime DG basis
% written by Pierson Guthrey

Nv1 = data.Nv1;
Nv2 = data.Nv2;
Nv3 = data.Nv3;
locs = data.locs;
v1centers = data.v1centers;
v2centers = data.v2centers;
v3centers = data.v3centers;
deltav1 = data.deltav1;
deltav2 = data.deltav2;
deltav3 = data.deltav3;

vectphi_Trans = data.vectphi_Trans;

DGcoeffs = zeros(data.thetaT,data.Nv1,data.Nv2,data.Nv3);

temp = zeros(data.thetaT,Nv1,Nv2,Nv3);


for iv1 = 1:Nv1
for iv2 = 1:Nv2
for iv3 = 1:Nv3
    cellcenter = [v1centers(iv1) v2centers(iv2) v3centers(iv3)];
    for k = 1:data.Pdp1
        tquad = data.Dp1list(k,1); 
        %tloc = data.Dp1quadlocs(k,1);
        v1quad = data.Dp1list(k,2); 
        v1loc = data.Dp1quadlocs(k,2);
        v1center = cellcenter(1);
        if data.space_dims >= 2
            v2quad = data.Dp1list(k,3);
            v2loc = data.Dp1quadlocs(k,3);
            v2center = cellcenter(2);
        else
            v2quad = 1;
            v2loc = inf;
        end
        if data.space_dims >= 3
            v3quad = data.Dp1list(k,4);
            v3loc = data.Dp1quadlocs(k,4);
            v3center = cellcenter(3);
        else
            v3quad = 1;
            v3loc = inf;
        end
        weight = data.Dp1quadwgts(k);        
        psik = data.vectpsi_Trans(:,:,tquad,v1quad,v2quad,v3quad);
        tau = locs(tquad);
        xii = locs(v1quad);
        eta = locs(v2quad);
        zta = locs(v3quad);
        v1tilde = v1centers(iv1)+deltav1*0.5*xii;
        v2tilde = v2centers(iv2)+deltav2*0.5*eta;
        v3tilde = v3centers(iv3)+deltav3*0.5*zta;
        ttilde = tnow+deltat/2*(tau+1);
        quadpoint = [ttilde,v1tilde,v2tilde,v3tilde]; 
        ftilde = f(quadpoint);
        temp(:,iv1,iv2,iv3) = temp(:,iv1,iv2,iv3) + psik*ftilde.*weight;
    end    
    DGcoeffs(:,iv1,iv2,iv3) = data.PSI_norm_inv*temp(:,iv1,iv2,iv3);
end
end
end

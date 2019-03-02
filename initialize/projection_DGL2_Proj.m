function [ DGcoeffs ] = projection_DGL2_Proj(f,data)
% Projects f onto the space DG basis
% written by Pierson Guthrey

theta = data.theta;
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

DGcoeffs = zeros(data.theta,data.Nv1,data.Nv2,data.Nv3);

temp = zeros(theta,Nv1,Nv2,Nv3);


data.PHI_norm_inv = 1/2^data.space_dims;
for iv1 = 1:Nv1
for iv2 = 1:Nv2
for iv3 = 1:Nv3
    for k = 1:data.Pd
        weight = data.Dquadwgts(k);
        v1quad = data.Dlist(k,1);
        if data.space_dims >= 2
            v2quad = data.Dlist(k,2);
        else
            v2quad = 1;
        end
        if data.space_dims >= 3
            v3quad = data.Dlist(k,3);
        else
            v3quad = 1;
        end
        phik = vectphi_Trans(:,:,v1quad,v2quad,v3quad);
        xii = locs(v1quad);
        eta = locs(v2quad);
        zta = locs(v3quad);
        v1tilde = v1centers(iv1)+deltav1*0.5*xii;
        v2tilde = v2centers(iv2)+deltav2*0.5*eta;
        v3tilde = v3centers(iv3)+deltav3*0.5*zta;
        quadpoint = [v1tilde,v2tilde,v3tilde]; 
        ftilde = f(quadpoint);
        temp(:,iv1,iv2,iv3) = temp(:,iv1,iv2,iv3) + phik*ftilde.*weight;
    end    
    DGcoeffs(:,iv1,iv2,iv3) = data.PHI_norm_inv*temp(:,iv1,iv2,iv3);
end
end
end

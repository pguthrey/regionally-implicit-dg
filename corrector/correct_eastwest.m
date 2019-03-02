function [correct_cell,maxspeedF,data] = correct_eastwest(qstar_1D,qeast_1D,qwest_1D,auxiliary,data,cellcenter)

correct_cell = 0;
v1quad = 1;
v1loc = inf;
v2quad = 1;
v2loc = inf;

error('whoops')

data.fluxeast(:,iv1) = zeros(data.theta,1);
data.fluxwest(:,iv1) = zeros(data.theta,1);

for k = 1:data.Pd
    weight = data.Dquadwgts(k);
    tquad = data.Dlist(k,1);
    if data.space_dims >= 2
        v1quad = data.Dlist(k,2);
        v1loc = data.Dquadlocs(k,1);    
        if data.space_dims >= 3
            v2quad = data.Dlist(k,3);
            v2loc = data.Dquadlocs(k,2);    
        end
    end
    %Flux in the x direction
    %Eastern Flux
    wstar = data.vectpsi_east(:,:,tquad,v1quad,v2quad)*qstar_1D;
    weast = data.vectpsi_west(:,:,tquad,v1quad,v2quad)*qeast_1D;
    phik = data.vectphi_east_Trans(:,:,v1quad,v2quad);
    v1 = cellcenter(1)+data.deltav1/2;
    v2 = cellcenter(2)+data.deltav2/2*v1loc;
    v3 = cellcenter(3)+data.deltav3/2*v2loc;
    quadpoint = [v1 v2 v3];
    [ Flux,speedeast ] = compute_F_numerical_flux(wstar,weast,quadpoint,data.appdata);
    correct_cell = correct_cell - (data.nuv1*weight)*phik*Flux;
    data.fluxeast(:,iv1) = data.fluxeast(:,iv1) - (data.nuv1*weight)*phik*Flux;    
    %Western Flux
    wstar = data.vectpsi_west(:,:,tquad,v1quad,v2quad)*qstar_1D;
    wwest = data.vectpsi_east(:,:,tquad,v1quad,v2quad)*qwest_1D;
    v1 = cellcenter(1)-data.deltav1/2;
    v2 = cellcenter(2)+data.deltav2/2*v1loc;
    v3 = cellcenter(3)+data.deltav3/2*v2loc;
    quadpoint = [v1 v2 v3];
    phik = data.vectphi_west_Trans(:,:,v1quad,v2quad);
    [ Flux, speedwest ] = compute_F_numerical_flux(wwest,wstar,quadpoint,data.appdata);
    correct_cell = correct_cell + (data.nuv1*weight).*phik*Flux;
    maxspeedF = max([speedeast speedwest]);
    
    data.fluxwest(:,iv1) = data.fluxwest(:,iv1) + (data.nuv1*weight).*phik*Flux;
    
end
for k = 1:data.Pdp1
    tquad = data.Dp1list(k,1);
    v1quad = data.Dp1list(k,2);
    v1loc = data.Dp1quadlocs(k,1); 
    if data.space_dims >= 2
        v2quad = data.Dp1list(k,3);
        v2loc = data.Dp1quadlocs(k,1); 
    else
        v2quad = 1;
        v2loc = inf;
    end
    if data.space_dims >= 3
        v3quad = data.Dp1list(k,4);
        v3loc = data.Dp1quadlocs(k,1); 
    else
        v3quad = 1;
        v3loc = inf;
    end
    weight = data.Dp1quadwgts(k);
    v1 = cellcenter(1)+data.deltav1/2*v1loc;
    v2 = cellcenter(2)+data.deltav2/2*v2loc;
    v3 = cellcenter(3)+data.deltav3/2*v3loc;
    quadpoint = [v1 v2 v3];
    wstar = data.vectpsi(:,:,tquad,v1quad,v2quad,v3quad)*qstar_1D;
    phikxii = data.vectphi_dxii_Trans(:,:,v1quad,v2quad,v3quad);
    fterm = phikxii*problem_F(wstar,quadpoint,data.appdata);
    correct_cell = correct_cell + data.nuv1*fterm*weight;
end

keyboard 

end


function [residual,maxspeed] = predictor_flux_east(qstar,qeast,data,cellcenter)
% written by Pierson Guthrey

residual = zeros(data.theta,1);

maxspeed = 0;

v1center = cellcenter(1);
v2center = cellcenter(2);
v3center = cellcenter(3);

for k = 1:data.Pdm1
    if data.space_dims >= 2
        v1quad = data.Dm1list(k,1);
        v1loc = data.Dm1quadlocs(k,1);
    else
        v1quad = 1;
        v1loc = inf;
    end
    if data.space_dims >= 3
        v2quad = data.Dm1list(k,2);
        v2loc = data.Dm1quadlocs(k,2);
    else
        v2quad = 1;
        v2loc = inf;
    end    
    weight = data.Dm1quadwgts(k);      

    wstar = data.vectphi_east(:,:,v1quad,v2quad)*qstar;
    weast = data.vectphi_west(:,:,v1quad,v2quad)*qeast;

    %East Boundary terms
    v1 = v1center-data.deltav1/2;        
    v2 = v2center+data.deltav2/2*v1loc;         
    v3 = v3center+data.deltav3/2*v2loc;   
    quadpoint = [v1 v2 v3];
    phik = data.vectphi_east_Trans(:,:,v1quad,v2quad);
    [ Flux,speed] = predictor_F_numerical_Flux(wstar,weast,quadpoint,data.appdata);
    residual = residual + data.nuv1*(phik*Flux)*weight;  

    maxspeed = max([maxspeed speed]);
end

end
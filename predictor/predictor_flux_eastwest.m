function [flux,maxspeedF] = predictor_flux_eastwest(qstar,data,cellcenter)
% written by Pierson Guthrey


error('huh?')

flux = 0; 
maxspeedF = 0;

v1center = cellcenter(1);
v2center = cellcenter(2);
v3center = cellcenter(3);

for k = 1:data.Pd
    tquad = data.Dlist(k,1); 
    %tloc = data.Dquadlocs(k,1);
    if data.space_dims >= 2
        v1quad = data.Dlist(k,2);
        v1loc = data.Dquadlocs(k,2);
    else
        v1quad = 1;
        v1loc = inf;
    end
    if data.space_dims >= 3
        v2quad = data.Dlist(k,3);
        v2loc = data.Dquadlocs(k,3);
    else
        v2quad = 1;
        v2loc = inf;
    end    
    weight = data.Dquadwgts(k);    
    wwest = data.vectpsi_east(:,:,tquad,v1quad,v2quad)*qstar;
    weast = data.vectpsi_west(:,:,tquad,v1quad,v2quad)*qeast;

    %East Boundary terms
    v1 = v1center+data.deltav1/2;        
    v2 = v2center+data.deltav2/2*v1loc;         
    v3 = v3center+data.deltav3/2*v2loc;   
    quadpoint = [v1 v2 v3];
    wstar = data.vectpsi_east(:,:,tquad,v1quad,v2quad)*qstar;
    psik = data.vectpsi_east_Trans(:,:,tquad,v1quad,v2quad);
    [ Flux,maxspeedF] = predictor_F_numerical_Flux(wstar,weast,quadpoint,data.appdata);
    fluxeast_temp = fluxeast_temp + data.nuv1*(psik*Flux)*weight;  

end


end
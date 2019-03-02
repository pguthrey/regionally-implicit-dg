function [correct_cell,maxspeedH] = correct_upprdown(correct_cell,qstar_1D,quppr_1D,qdown_1D,auxiliary,data,cellcenter)


error('whoops')


for k = 1:data.Pd
    weight = data.Dquadwgts(k);
    tquad = data.Dlist(k,1);
    if data.space_dims >= 2
        v1quad = data.Dlist(k,2);
        v1loc = data.Dquadlocs(k,1);    
    else
        v1quad = 1;
        v1loc = inf;
    end
    if data.space_dims >= 3
        v2quad = data.Dlist(k,3);
        v2loc = data.Dquadlocs(k,2);    
    else
        v2quad = 1;
        v2loc = inf;
    end
    %Flux in the z direction
    %Up Flux
    wstar = data.vectpsi_uppr(:,:,tquad,v1quad,v2quad)*qstar_1D;
    wuppr = data.vectpsi_down(:,:,tquad,v1quad,v2quad)*quppr_1D;
    v1 = cellcenter(1)+data.deltav1/2*v1loc;
    v2 = cellcenter(2)+data.deltav2/2*v2loc;
    v3 = cellcenter(3)+data.deltav3/2;
    quadpoint = [v1 v2 v3];
    phik = data.vectphi_uppr_Trans(:,:,v1quad,v2quad);
    [ Flux,speeduppr ] = compute_H_numerical_flux(wstar,wuppr,quadpoint,data.appdata);
    correct_cell = correct_cell - (data.nuv3*weight)*phik*Flux;

    %Down Flux
    wstar = data.vectpsi_down(:,:,tquad,v1quad,v2quad)*qstar_1D;
    wdown = data.vectpsi_uppr(:,:,tquad,v1quad,v2quad)*qdown_1D;
    v1 = cellcenter(1)+data.deltav1/2*data.locs(v1quad);
    v2 = cellcenter(2)+data.deltav2/2*data.locs(v2quad);
    v3 = cellcenter(3)-data.deltav3/2;
    quadpoint = [v1 v2 v3];
    phik = data.vectphi_down_Trans(:,:,v1quad,v2quad);
    [ Flux, speeddown ] = compute_H_numerical_flux(wdown,wstar,quadpoint,data.appdata);
    correct_cell = correct_cell + (data.nuv3*weight)*phik*Flux;            

    maxspeedH = max([speeduppr speeddown]);       
end
for k = 1:data.Pdp1
    tquad = data.Dp1list(k,1);
    v1quad = data.Dp1list(k,2);
    v1loc = data.Dp1quadlocs(k,1); 
    v2quad = data.Dp1list(k,3);
    v2loc = data.Dp1quadlocs(k,1); 
    v3quad = data.Dp1list(k,4);
    v3loc = data.Dp1quadlocs(k,1); 
    weight = data.Dp1quadwgts(k);
    v1 = cellcenter(1)+data.deltav1/2*v1loc;
    v2 = cellcenter(2)+data.deltav2/2*v2loc;
    v3 = cellcenter(2)+data.deltav3/2*v3loc;
    quadpoint = [v1 v2 v3];
    wstar = data.vectpsi(:,:,tquad,v1quad,v2quad,v3quad)*qstar_1D;
    phikzta = data.vectphi_dzta_Trans(:,:,v1quad,v2quad,v3quad);
    hterm = phikzta*problem_H(wstar,quadpoint,data.appdata);
    correct_cell = correct_cell + data.nuv3*hterm*weight;
end
end


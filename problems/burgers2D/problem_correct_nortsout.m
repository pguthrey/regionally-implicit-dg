function [correct_cell,maxspeedG,data] = problem_correct_nortsout(correct_cell,qstar,qnort,qsout,data,cellcenter)

for k = 1:data.Pd
    weight = data.Dquadwgts(k);
    tquad = data.Dlist(k,1);
    v1quad = data.Dlist(k,2);
    v1loc = data.Dquadlocs(k,1);    
    if data.space_dims >= 3
        v2quad = data.Dlist(k,3);
        v2loc = data.Dquadlocs(k,2);    
    else
        v2quad = 1;
        v2loc = inf;
    end
    %Flux in the y direction
    %Northern Flux
    wstar = data.vectpsi_nort(:,:,tquad,v1quad,v2quad)*qstar;
    wnort = data.vectpsi_sout(:,:,tquad,v1quad,v2quad)*qnort;
    v1 = cellcenter(1)+data.deltav1/2*v1loc;
    v2 = cellcenter(2)+data.deltav2/2;
    v3 = cellcenter(3)+data.deltav3/2*v2loc;
    quadpoint = [v1 v2 v3];
    phik = data.vectphi_nort_Trans(:,:,v1quad,v2quad);
    [ Flux,speednort ] = compute_G_numerical_flux(wstar,wnort,quadpoint,data.appdata);
    correct_cell = correct_cell - (data.nuv2*weight)*phik*Flux;
    
    %Southern Flux
    wstar = data.vectpsi_sout(:,:,tquad,v1quad,v2quad)*qstar;
    wsout = data.vectpsi_nort(:,:,tquad,v1quad,v2quad)*qsout;
    v1 = cellcenter(1)+data.deltav1/2*v1loc;
    v2 = cellcenter(2)-data.deltav2/2;
    v3 = cellcenter(3)+data.deltav3/2*v2loc;
    quadpoint = [v1 v2 v3];
    phik = data.vectphi_sout_Trans(:,:,v1quad,v2quad);
    [ Flux, speedsout ] = compute_G_numerical_flux(wsout,wstar,quadpoint,data.appdata);
    correct_cell = correct_cell + (data.nuv2*weight)*phik*Flux;            
            
    maxspeedG = max([speednort speedsout]);

end
for k = 1:data.Pdp1
        tquad = data.Dp1list(k,1);
        v1quad = data.Dp1list(k,2);
        v1loc = data.Dp1quadlocs(k,1); 
        v2quad = data.Dp1list(k,3);
        v2loc = data.Dp1quadlocs(k,1); 
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
        v3 = cellcenter(2)+data.deltav3/2*v3loc;
        quadpoint = [v1 v2 v3];
        wstar = data.vectpsi(:,:,tquad,v1quad,v2quad,v3quad)*qstar;
        phiketa = data.vectphi_deta_Trans(:,:,v1quad,v2quad,v3quad);
        gterm = phiketa*problem_G(wstar,quadpoint,data.appdata);
        correct_cell = correct_cell + data.nuv2*gterm*weight;
end


end


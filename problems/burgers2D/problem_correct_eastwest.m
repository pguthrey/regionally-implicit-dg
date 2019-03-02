function [correct_cell,maxspeedF,data] = problem_correct_eastwest(qstar_1D,qeast_1D,qwest_1D,data,cellcenter,iv1)

correct_cell = 0;
v2quad = 1;
quadpoint = nan;

%data.fluxeast(:,iv1) = zeros(data.theta,1);
%data.fluxwest(:,iv1) = zeros(data.theta,1);

for k = 1:data.Pd
    weight = data.Dquadwgts(k);
    tquad = data.Dlist(k,1);
    v1quad = data.Dlist(k,2);
    %Flux in the x direction
    %Eastern Flux
    wstar = data.vectpsi_east(:,:,tquad,v1quad,v2quad)*qstar_1D;
    weast = data.vectpsi_west(:,:,tquad,v1quad,v2quad)*qeast_1D;
    phik = data.vectphi_east_Trans(:,:,v1quad,v2quad);
    [ Flux,speedeast ] = compute_F_numerical_flux(wstar,weast,quadpoint,data.appdata);
    correct_cell = correct_cell - (data.nuv1*weight)*phik*Flux;
    
    %Western Flux
    wstar = data.vectpsi_west(:,:,tquad,v1quad,v2quad)*qstar_1D;
    wwest = data.vectpsi_east(:,:,tquad,v1quad,v2quad)*qwest_1D;
    phik = data.vectphi_west_Trans(:,:,v1quad,v2quad);
    [ Flux, speedwest ] = compute_F_numerical_flux(wwest,wstar,quadpoint,data.appdata);
    correct_cell = correct_cell + (data.nuv1*weight).*phik*Flux;
    maxspeedF = max([speedeast speedwest]);
end
for k = 1:data.Pdp1
    tquad = data.Dp1list(k,1);
    v1quad = data.Dp1list(k,2);
    v2quad = data.Dp1list(k,3);
    weight = data.Dp1quadwgts(k);  
    wstar = data.vectpsi(:,:,tquad,v1quad,v2quad)*qstar_1D;
    phikxii = data.vectphi_dxii_Trans(:,:,v1quad,v2quad);
    fterm = phikxii*problem_F(wstar,quadpoint,data.appdata);
    correct_cell = correct_cell + data.nuv1*fterm*weight;
end
end

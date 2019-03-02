function [correct_cell,maxspeedF,data] = problem_correct_eastwest(qstar,qeast,qwest,data,cellcenter,iv1)

maxspeedF = 0;
correct_cell = 0;

vectphi_east_Trans = data.vectphi_east_Trans;
vectphi_west_Trans = data.vectphi_west_Trans;
vectphi_dxii_Trans = data.vectphi_dxii_Trans;

vectpsi_east = data.vectpsi_east;
vectpsi_west = data.vectpsi_west;
vectpsi = data.vectpsi;

Pd = data.Pd;
Pdp1 = data.Pdp1;
Dlist = data.Dlist;
space_dims = data.space_dims;
Dp1list = data.Dp1list;
Dp1quadwgts = data.Dp1quadwgts;
Dquadwgts = data.Dquadwgts;

nuv1 = data.nuv1;

for k = 1:Pd
    weight = Dquadwgts(k);
    tquad = Dlist(k,1);
    v1quad = Dlist(k,2);
    if space_dims >= 3
        v2quad = Dlist(k,3);
    else
        v2quad = 1;
    end
    %Flux in the y direction
    %Northern Flux
    
    wwest_p = vectpsi_east(:,:,tquad,v1quad,v2quad)*qwest;
    weast_m = vectpsi_west(:,:,tquad,v1quad,v2quad)*qeast;
    wstar_p = vectpsi_east(:,:,tquad,v1quad,v2quad)*qstar;
    wstar_m = vectpsi_west(:,:,tquad,v1quad,v2quad)*qstar;

    Fwest_p = data.appdata.nuf*wwest_p;
    Feast_m = data.appdata.nuf*weast_m;

    Fstar_p = data.appdata.nuf*wstar_p;    
    Fstar_m = data.appdata.nuf*wstar_m;
    
    speedeast = max(abs(data.appdata.nuf),abs(data.appdata.nuf));
    speedwest = max(abs(data.appdata.nuf),abs(data.appdata.nuf));
    
    %North Boundary terms
    psik = vectphi_east_Trans(:,:,v1quad,v2quad);
    Flux = (Fstar_p + Feast_m - speedeast*(weast_m-wstar_p))/2;%  compute_F_numerical_flux(wstar,weast,quadpoint,appdata);
    correct_cell = correct_cell + (psik*Flux)*(nuv1*weight);  
    
    %west boundary terms

    psik = vectphi_west_Trans(:,:,v1quad,v2quad);
    Flux =  (Fstar_m + Fwest_p - speedwest*(wstar_m-wwest_p))/2;%compute_F_numerical_flux(wwest,wstar,quadpoint,appdata);
    correct_cell = correct_cell - (psik*Flux)*(nuv1*weight);
    
    maxspeedF_quad = max(speedeast,speedwest);
    maxspeedF = max(maxspeedF,maxspeedF_quad);    
    
end
for k = 1:Pdp1
        tquad = Dp1list(k,1);
        v1quad = Dp1list(k,2);
        v2quad = Dp1list(k,3);
        if space_dims >= 3
            v3quad = Dp1list(k,4);
        else
            v3quad = 1;
        end
        weight = Dp1quadwgts(k);
        wstar = vectpsi(:,:,tquad,v1quad,v2quad,v3quad)*qstar;
        Fstar = data.appdata.nuf*wstar;
        phikxii = vectphi_dxii_Trans(:,:,v1quad,v2quad,v3quad);
        correct_cell = correct_cell - nuv1*weight*phikxii*Fstar;
end

%disp(['correct eastwest speed = ' num2str(maxspeedF)])


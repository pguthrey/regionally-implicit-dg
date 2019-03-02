function [correct_cell,auxiliary,maxspeedG,data] = problem_correct_nortsout(correct_cell,qstar,qnort,qsout,auxiliary,data,cellcenter,cellindex)

vectphi_nort_Trans = data.vectphi_nort_Trans;
vectphi_sout_Trans = data.vectphi_sout_Trans;
vectphi_deta_Trans = data.vectphi_deta_Trans;

vectpsi_nort = data.vectpsi_nort;
vectpsi_sout = data.vectpsi_sout;
vectpsi = data.vectpsi;

Pd = data.Pd;
Pdp1 = data.Pdp1;
Dlist = data.Dlist;
space_dims = data.space_dims;
Dp1list = data.Dp1list;
Dp1quadwgts = data.Dp1quadwgts;
Dquadwgts = data.Dquadwgts;

nuv2 = data.nuv2;
maxspeedG = 0;

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
    
    wsout_p = vectpsi_nort(:,:,tquad,v1quad,v2quad)*qsout;
    wnort_m = vectpsi_sout(:,:,tquad,v1quad,v2quad)*qnort;
    wstar_p = vectpsi_nort(:,:,tquad,v1quad,v2quad)*qstar;
    wstar_m = vectpsi_sout(:,:,tquad,v1quad,v2quad)*qstar;

    Fsout_p = wsout_p^2/2;
    Fnort_m = wnort_m^2/2;

    Fstar_p = wstar_p^2/2;    
    Fstar_m = wstar_m^2/2;
    
    speednort = max(abs(wstar_p),abs(wnort_m));
    speedsout = max(abs(wstar_m),abs(wsout_p));
    
    %North Boundary terms
    psik = vectphi_nort_Trans(:,:,v1quad,v2quad);
    Flux = (Fstar_p + Fnort_m - speednort*(wnort_m-wstar_p))/2;%  compute_G_numerical_flux(wstar,wnort,quadpoint,appdata);
    correct_cell = correct_cell + (psik*Flux)*(nuv2*weight);  
    
    %sout boundary terms

    psik = vectphi_sout_Trans(:,:,v1quad,v2quad);
    Flux =  (Fstar_m + Fsout_p - speedsout*(wstar_m-wsout_p))/2;%compute_G_numerical_flux(wsout,wstar,quadpoint,appdata);
    correct_cell = correct_cell - (psik*Flux)*(nuv2*weight);
    
    maxspeedG_quad = max(speednort,speedsout);
    maxspeedG = max(maxspeedG,maxspeedG_quad);    
    
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
        phiketa = vectphi_deta_Trans(:,:,v1quad,v2quad,v3quad);
        gterm = phiketa*wstar^2/2;
        correct_cell = correct_cell - nuv2*gterm*weight;
end


end


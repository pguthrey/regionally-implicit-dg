function [correct_cell,auxiliary,maxspeedF,data] = problem_correct_eastwest(qstar,qeast,qwest,auxiliary,data,cellcenter,cellindex)

vectphi_east = data.vectphi_east;
vectphi_west = data.vectphi_west;

vectphi_east_Trans = data.vectphi_east_Trans;
vectphi_west_Trans = data.vectphi_west_Trans;
vectphi_dxii_Trans = data.vectphi_dxii_Trans;

vectpsi_east = data.vectpsi_east;
vectpsi_west = data.vectpsi_west;
vectphi = data.vectphi;
vectpsi = data.vectpsi;

Pd = data.Pd;
Pdp1 = data.Pdp1;
Dlist = data.Dlist;
space_dims = data.space_dims;
Dp1list = data.Dp1list;
Dp1quadwgts = data.Dp1quadwgts;
Dquadwgts = data.Dquadwgts;

nuv1 = data.nuv1;
maxspeedF = 0;

ivx = cellindex(1);
ivy = cellindex(2);
ivz = cellindex(3);
periodic = @(i,N) mod(i-1,N)+1; 
Nv1 = data.Nv1;
Nv2 = data.Nv2;
Nv3 = data.Nv3;
Pwest = auxiliary.total_pressure(:,periodic(ivx-1,Nv1),ivy,ivz);
Pstar = auxiliary.total_pressure(:,periodic(ivx  ,Nv1),ivy,ivz);
Peast = auxiliary.total_pressure(:,periodic(ivx+1,Nv1),ivy,ivz);

correct_cell = 0;

Flux_x =@(q,p)[ q(2) ; 
              q(2)/q(1)*q(2) + p ;
              q(2)/q(1)*q(3)  ;
              q(2)/q(1)*q(4)  ;
              q(2)/q(1)*q(5) ;
              q(2)/q(1)*q(6) ;
              q(2)/q(1)*q(7) ]; 

DPDNeast = auxiliary.dpdn(:,periodic(ivx+1,Nv1),periodic(ivy,Nv2),periodic(ivz,Nv3));
DPDNstar = auxiliary.dpdn(:,periodic(ivx ,Nv1),periodic(ivy,Nv2),periodic(ivz,Nv3));
DPDNwest = auxiliary.dpdn(:,periodic(ivx-1,Nv1),periodic(ivy,Nv2),periodic(ivz,Nv3));

speed_x = @(q,dpdn) abs(q(2)/q(1)) + abs(sqrt(dpdn)) ; 

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

    pwest_p = vectphi_east(1,1:data.Ls,v1quad,v2quad)*Pwest;
    peast_m = vectphi_west(1,1:data.Ls,v1quad,v2quad)*Peast;
    pstar_p = vectphi_east(1,1:data.Ls,v1quad,v2quad)*Pstar;
    pstar_m = vectphi_west(1,1:data.Ls,v1quad,v2quad)*Pstar;

    dpdnwest_p = vectphi_east(1,1:data.Ls,v1quad,v2quad)*DPDNwest;
    dpdneast_m = vectphi_west(1,1:data.Ls,v1quad,v2quad)*DPDNeast;
    dpdnstar_p = vectphi_east(1,1:data.Ls,v1quad,v2quad)*DPDNstar;
    dpdnstar_m = vectphi_west(1,1:data.Ls,v1quad,v2quad)*DPDNstar;
    
    Fwest_p = Flux_x(wwest_p,pwest_p);
    Feast_m = Flux_x(weast_m,peast_m);

    Fstar_p = Flux_x(wstar_p,pstar_p);    
    Fstar_m = Flux_x(wstar_m,pstar_m);
        
    %North Boundary terms
    
    speedleft=speed_x(wstar_p, dpdnstar_p);
    speedavg=speed_x((wstar_p+weast_m)/2, (dpdnstar_p+dpdneast_m)/2 ) ;
    speedrite=speed_x(weast_m, dpdneast_m);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);    
    psik = vectphi_east_Trans(:,:,v1quad,v2quad);
    Flux = (Fstar_p + Feast_m - speed*(weast_m-wstar_p))/2;%  compute_G_numerical_flux(wstar,weast,quadpoint,appdata);
    correct_cell = correct_cell + (psik*Flux)*(nuv1*weight);  
    
    %west boundary terms
    speedleft=speed_x(wwest_p, dpdnwest_p);
    speedavg=speed_x((wwest_p+wstar_m)/2, (dpdnwest_p+dpdnstar_m)/2 ) ;
    speedrite=speed_x(wstar_m, dpdnstar_m);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);    

    psik = vectphi_west_Trans(:,:,v1quad,v2quad);
    Flux =  (Fstar_m + Fwest_p - speed*(wstar_m-wwest_p))/2;%compute_G_numerical_flux(wwest,wstar,quadpoint,appdata);
    correct_cell = correct_cell - (psik*Flux)*(nuv1*weight);
    
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
        pstar = vectphi(1,1:data.Ls,v1quad,v2quad,v3quad)*Pstar;
        phikxii = vectphi_dxii_Trans(:,:,v1quad,v2quad,v3quad);
        gterm = phikxii*Flux_x(wstar,pstar);
        correct_cell = correct_cell - nuv1*gterm*weight;
end


end


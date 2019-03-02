function [correct_cell,auxiliary,maxspeedF,data] = problem_correct_nortsout(correct_cell,qstar,qnort,qsout,auxiliary,data,cellcenter,cellindex)

vectphi_nort = data.vectphi_nort;
vectphi_sout = data.vectphi_sout;

vectphi_nort_Trans = data.vectphi_nort_Trans;
vectphi_sout_Trans = data.vectphi_sout_Trans;
vectphi_deta_Trans = data.vectphi_deta_Trans;

vectpsi_nort = data.vectpsi_nort;
vectpsi_sout = data.vectpsi_sout;
vectphi = data.vectphi;
vectpsi = data.vectpsi;

Pd = data.Pd;
Pdp1 = data.Pdp1;
Dlist = data.Dlist;
space_dims = data.space_dims;
Dp1list = data.Dp1list;
Dp1quadwgts = data.Dp1quadwgts;
Dquadwgts = data.Dquadwgts;

nuv2 = data.nuv2;
maxspeedF = 0;

ivx = cellindex(1);
ivy = cellindex(2);
ivz = cellindex(3);
periodic = @(i,N) mod(i-1,N)+1; 
Nv1 = data.Nv1;
Nv2 = data.Nv2;
Nv3 = data.Nv3;
Psout = auxiliary.total_pressure(:,ivx,periodic(ivy-1,Nv2),ivz);
Pstar = auxiliary.total_pressure(:,ivx,periodic(ivy  ,Nv2),ivz);
Pnort = auxiliary.total_pressure(:,ivx,periodic(ivy+1,Nv2),ivz);

DPDNnort = auxiliary.dpdn(:,periodic(ivx,Nv1),periodic(ivy+1,Nv2),periodic(ivz,Nv3));
DPDNstar = auxiliary.dpdn(:,periodic(ivx,Nv1),periodic(ivy,Nv2),periodic(ivz,Nv3));
DPDNsout = auxiliary.dpdn(:,periodic(ivx,Nv1),periodic(ivy-1,Nv2),periodic(ivz,Nv3));


Flux_y = @(q,p)[ q(3);
 (q(2)*q(3))/q(1);
    q(3)^2/q(1) + p ;
 (q(3)*q(4))/q(1);
 (q(5)*q(3))/q(1);
 (q(6)*q(3))/q(1);
 (q(7)*q(3))/q(1);];

speed_y = @(q,dpdn) abs(q(3)/q(1)) + abs(sqrt(dpdn)) ; 

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

    psout_p = vectphi_nort(1,1:data.Ls,v1quad,v2quad)*Psout;
    pnort_m = vectphi_sout(1,1:data.Ls,v1quad,v2quad)*Pnort;
    pstar_p = vectphi_nort(1,1:data.Ls,v1quad,v2quad)*Pstar;
    pstar_m = vectphi_sout(1,1:data.Ls,v1quad,v2quad)*Pstar;

    dpdnsout_p = vectphi_nort(1,1:data.Ls,v1quad,v2quad)*DPDNsout;
    dpdnnort_m = vectphi_sout(1,1:data.Ls,v1quad,v2quad)*DPDNnort;
    dpdnstar_p = vectphi_nort(1,1:data.Ls,v1quad,v2quad)*DPDNstar;
    dpdnstar_m = vectphi_sout(1,1:data.Ls,v1quad,v2quad)*DPDNstar;
    
    Fsout_p = Flux_y(wsout_p,psout_p);
    Fnort_m = Flux_y(wnort_m,pnort_m);

    Fstar_p = Flux_y(wstar_p,pstar_p);    
    Fstar_m = Flux_y(wstar_m,pstar_m);
    
    speednort = max(abs(wstar_p),abs(wnort_m));
    speedsout = max(abs(wstar_m),abs(wsout_p));
    
    %North Boundary terms
    
    speedleft=speed_y(wstar_p, dpdnstar_p);
    speedavg=speed_y((wstar_p+wnort_m)/2, (dpdnstar_p+dpdnnort_m)/2 ) ;
    speedrite=speed_y(wnort_m, dpdnnort_m);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);    
    psik = vectphi_nort_Trans(:,:,v1quad,v2quad);
    Flux = (Fstar_p + Fnort_m - speed*(wnort_m-wstar_p))/2;%  compute_G_numerical_flux(wstar,wnort,quadpoint,appdata);
    correct_cell = correct_cell + (psik*Flux)*(nuv2*weight);  
    
    %sout boundary terms
    speedleft=speed_y(wsout_p, dpdnsout_p);
    speedavg=speed_y((wsout_p+wstar_m)/2, (dpdnsout_p+dpdnstar_m)/2 ) ;
    speedrite=speed_y(wstar_m, dpdnstar_m);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);    

    psik = vectphi_sout_Trans(:,:,v1quad,v2quad);
    Flux =  (Fstar_m + Fsout_p - speed*(wstar_m-wsout_p))/2;%compute_G_numerical_flux(wsout,wstar,quadpoint,appdata);
    correct_cell = correct_cell - (psik*Flux)*(nuv2*weight);
    
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
        phiketa = vectphi_deta_Trans(:,:,v1quad,v2quad,v3quad);
        gterm = phiketa*Flux_y(wstar,pstar);
        correct_cell = correct_cell - nuv2*gterm*weight;
end



end


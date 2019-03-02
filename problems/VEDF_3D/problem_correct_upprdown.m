function [correct_cell,auxiliary,maxspeedF,data] = problem_correct_upprdown(correct_cell,qstar,quppr,qdown,auxiliary,data,cellcenter,cellindex)

vectphi_uppr = data.vectphi_uppr;
vectphi_down = data.vectphi_down;

vectphi_uppr_Trans = data.vectphi_uppr_Trans;
vectphi_down_Trans = data.vectphi_down_Trans;
vectphi_deta_Trans = data.vectphi_deta_Trans;

vectpsi_uppr = data.vectpsi_uppr;
vectpsi_down = data.vectpsi_down;
vectphi = data.vectphi;
vectpsi = data.vectpsi;

Pd = data.Pd;
Pdp1 = data.Pdp1;
Dlist = data.Dlist;
space_dims = data.space_dims;
Dp1list = data.Dp1list;
Dp1quadwgts = data.Dp1quadwgts;
Dquadwgts = data.Dquadwgts;

nuv3 = data.nuv3;
maxspeedF = 0;

ivx = cellindex(1);
ivy = cellindex(2);
ivz = cellindex(3);
periodic = @(i,N) mod(i-1,N)+1; 
Nv1 = data.Nv1;
Nv2 = data.Nv2;
Nv3 = data.Nv3;
Pdown = auxiliary.total_pressure(:,ivx,ivy,periodic(ivz-1,Nv3));
Pstar = auxiliary.total_pressure(:,ivx,ivy,periodic(ivz  ,Nv3));
Puppr = auxiliary.total_pressure(:,ivx,ivy,periodic(ivz+1,Nv3));

DPDNuppr = auxiliary.dpdn(:,periodic(ivx,Nv1),periodic(ivy,Nv2),periodic(ivz+1,Nv3));
DPDNstar= auxiliary.dpdn(:,periodic(ivx ,Nv1),periodic(ivy,Nv2),periodic(ivz,Nv3));
DPDNdown = auxiliary.dpdn(:,periodic(ivx,Nv1),periodic(ivy,Nv2),periodic(ivz-1,Nv3));

 
Flux_z = @(q,p)[ q(4);
 (q(2)*q(4))/q(1);
 (q(3)*q(4))/q(1);
    q(4)^2/q(1) + p ;
 (q(5)*q(4))/q(1);
 (q(6)*q(4))/q(1);
 (q(7)*q(4))/q(1);];

speed_z = @(q,dpdn) abs(q(4)/q(1)) + abs(sqrt(dpdn)) ; 


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
    wdown_p = vectpsi_uppr(:,:,tquad,v1quad,v2quad)*qdown;
    wuppr_m = vectpsi_down(:,:,tquad,v1quad,v2quad)*quppr;
    wstar_p = vectpsi_uppr(:,:,tquad,v1quad,v2quad)*qstar;
    wstar_m = vectpsi_down(:,:,tquad,v1quad,v2quad)*qstar;

    pdown_p = vectphi_uppr(1,1:data.Ls,v1quad,v2quad)*Pdown;
    puppr_m = vectphi_down(1,1:data.Ls,v1quad,v2quad)*Puppr;
    pstar_p = vectphi_uppr(1,1:data.Ls,v1quad,v2quad)*Pstar;
    pstar_m = vectphi_down(1,1:data.Ls,v1quad,v2quad)*Pstar;

    dpdndown_p = vectphi_uppr(1,1:data.Ls,v1quad,v2quad)*DPDNdown;
    dpdnuppr_m = vectphi_down(1,1:data.Ls,v1quad,v2quad)*DPDNuppr;
    dpdnstar_p = vectphi_uppr(1,1:data.Ls,v1quad,v2quad)*DPDNstar;
    dpdnstar_m = vectphi_down(1,1:data.Ls,v1quad,v2quad)*DPDNstar;
    
    Fdown_p = Flux_z(wdown_p,pdown_p);
    Fuppr_m = Flux_z(wuppr_m,puppr_m);

    Fstar_p = Flux_z(wstar_p,pstar_p);    
    Fstar_m = Flux_z(wstar_m,pstar_m);
  
    
    %North Boundary terms
    
    speedleft=speed_z(wstar_p, dpdnstar_p);
    speedavg=speed_z((wstar_p+wuppr_m)/2, (dpdnstar_p+dpdnuppr_m)/2 ) ;
    speedrite=speed_z(wuppr_m, dpdnuppr_m);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);    
    psik = vectphi_uppr_Trans(:,:,v1quad,v2quad);
    Flux = (Fstar_p + Fuppr_m - speed*(wuppr_m-wstar_p))/2;%  compute_G_numerical_flux(wstar,wuppr,quadpoint,appdata);
    correct_cell = correct_cell + (psik*Flux)*(nuv3*weight);  
    
    %down boundary terms
    speedleft=speed_z(wdown_p, dpdndown_p);
    speedavg=speed_z((wdown_p+wstar_m)/2, (dpdndown_p+dpdnstar_m)/2 ) ;
    speedrite=speed_z(wstar_m, dpdnstar_m);
    speed=max(speedleft,speedrite);
    speed=max(speed,speedavg);
    maxspeedF=max(maxspeedF,speed);    

    psik = vectphi_down_Trans(:,:,v1quad,v2quad);
    Flux =  (Fstar_m + Fdown_p - speed*(wstar_m-wdown_p))/2;%compute_G_numerical_flux(wdown,wstar,quadpoint,appdata);
    correct_cell = correct_cell - (psik*Flux)*(nuv3*weight);
        
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
        gterm = phiketa*Flux_z(wstar,pstar);
        correct_cell = correct_cell - nuv3*gterm*weight;
end


end


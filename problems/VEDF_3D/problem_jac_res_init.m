function [Jacobian,residual,auxiliary] ...
          = problem_jac_res_init(qstar,qpast,Pstar,dPstardn,total_viscocity,auxiliary,data)
% Forms the %Jacobian of equations given by the DG data.method for a fiv1ed cell
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    DGsubregion
%           prevcell
%           starcell
%           %Jac_truncdirs
% OUTdata.PUTS   %Jacobian
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------    

residual = (data.I_futr - data.Psi_dtau)*qstar-data.I_past*qpast;
Jacobian = (data.I_futr - data.Psi_dtau);

vectpsi_dxii_Trans = data.vectpsi_dxii_Trans;
vectpsi_deta_Trans = data.vectpsi_deta_Trans;
vectpsi_dzta_Trans = data.vectpsi_dzta_Trans;
vectpsi = data.vectpsi;
vectpsi_Trans = data.vectpsi_Trans;

v2quad = 1;
v3quad = 1;
Dp1list = data.Dp1list;
Pdp1 = data.Pdp1;
space_dims = data.space_dims;
Dp1quadwgts = data.Dp1quadwgts;
nuv1 = data.nuv1;
nuv2 = data.nuv2;
nuv3 = data.nuv3;

Flux_x = @(q,p)[ q(2);
    q(2)^2/q(1) + p;
 (q(2)*q(3))/q(1);
 (q(2)*q(4))/q(1);
 (q(5)*q(2))/q(1);
 (q(6)*q(2))/q(1);
 (q(7)*q(2))/q(1);];
 
Flux_y = @(q,p)[ q(3);
 (q(2)*q(3))/q(1);
    q(3)^2/q(1) + p ;
 (q(3)*q(4))/q(1);
 (q(5)*q(3))/q(1);
 (q(6)*q(3))/q(1);
 (q(7)*q(3))/q(1);];
 
Flux_z = @(q,p)[ q(4);
 (q(2)*q(4))/q(1);
 (q(3)*q(4))/q(1);
    q(4)^2/q(1) + p ;
 (q(5)*q(4))/q(1);
 (q(6)*q(4))/q(1);
 (q(7)*q(4))/q(1);];

tau = 1;

% EXTERNAL FORCES SET TO ZERO
Source = @(q)[  0;
                q(5)/q(1);
                q(6)/q(1);
                q(7)/q(1);
                (total_viscocity(1)-q(5))/tau
                (total_viscocity(2)-q(6))/tau
                (total_viscocity(3)-q(7))/tau ];
DSourcedQ = @(q)[  0,0,0,0,0,0,0;
                -q(5)/q(1)^2,0,0,0,1/q(1),0 ,0;
                -q(6)/q(1)^2,0,0,0,0,1/q(1) ,0;
                -q(7)/q(1)^2,0,0,0,0,0,1/q(1) ;
                0,0,0,0,1/tau,0 ,0;
                0,0,0,0,0,1/tau ,0;
                0,0,0,0,0,0,1/tau ;];
            
            
DFlux_dq_x = @(q,dpdn) [[               0,        1,    0,    0,    0,    0,    0];
[ dpdn - q(2)^2/q(1)^2, (2*q(2))/q(1),    0,    0,    0,    0,    0];
[    -(q(2)*q(3))/q(1)^2,     q(3)/q(1), q(2)/q(1),    0,    0,    0,    0];
[    -(q(2)*q(4))/q(1)^2,     q(4)/q(1),    0, q(2)/q(1),    0,    0,    0];
[    -(q(5)*q(2))/q(1)^2,     q(5)/q(1),    0,    0, q(2)/q(1),    0,    0];
[    -(q(6)*q(2))/q(1)^2,     q(6)/q(1),    0,    0,    0, q(2)/q(1),    0];
[    -(q(7)*q(2))/q(1)^2,     q(7)/q(1),    0,    0,    0,    0, q(2)/q(1)]];


DFlux_dq_y = @(q,dpdn) [[               0,    0,        1,    0,    0,    0,    0];
[    -(q(2)*q(3))/q(1)^2, q(3)/q(1),     q(2)/q(1),    0,    0,    0,    0];
[ dpdn - q(3)^2/q(1)^2,    0, (2*q(3))/q(1),    0,    0,    0,    0];
[    -(q(3)*q(4))/q(1)^2,    0,     q(4)/q(1), q(3)/q(1),    0,    0,    0];
[    -(q(5)*q(3))/q(1)^2,    0,     q(5)/q(1),    0, q(3)/q(1),    0,    0];
[    -(q(6)*q(3))/q(1)^2,    0,     q(6)/q(1),    0,    0, q(3)/q(1),    0];
[    -(q(7)*q(3))/q(1)^2,    0,     q(7)/q(1),    0,    0,    0, q(3)/q(1)]];

DFlux_dq_z = @(q,dpdn) [[               0,    0,    0,        1,    0,    0,    0];
[    -(q(2)*q(4))/q(1)^2, q(4)/q(1),    0,     q(2)/q(1),    0,    0,    0];
[    -(q(3)*q(4))/q(1)^2,    0, q(4)/q(1),     q(3)/q(1),    0,    0,    0];
[ dpdn - q(4)^2/q(1)^2,    0,    0, (2*q(4))/q(1),    0,    0,    0];
[    -(q(5)*q(4))/q(1)^2,    0,    0,     q(5)/q(1), q(4)/q(1),    0,    0];
[    -(q(6)*q(4))/q(1)^2,    0,    0,     q(6)/q(1),    0, q(4)/q(1),    0];
[    -(q(7)*q(4))/q(1)^2,    0,    0,     q(7)/q(1),    0,    0, q(4)/q(1)];];



for k = 1:Pdp1
        tquad = Dp1list(k,1); 
        %tloc = Dp1quadlocs(k,1);
        v1quad = Dp1list(k,2); 
        if space_dims >= 2
            v2quad = Dp1list(k,3);
            if space_dims >= 3
                v3quad = Dp1list(k,4);
            end
        end
        weight = Dp1quadwgts(k);
        %Derivative terms
        
        %eastwest
        wstar = vectpsi(:,:,tquad,v1quad,v2quad,v3quad)*qstar;
        pstar = data.vectphi(1,1:data.Ls,v1quad,v2quad,v3quad)*Pstar;
        dpstardn = data.vectphi(1,1:data.Ls,v1quad,v2quad,v3quad)*dPstardn;
        
        psikdxii = vectpsi_dxii_Trans(:,:,tquad,v1quad,v2quad,v3quad);
        Fstar = Flux_x(wstar,pstar);
        fterm = (nuv1*weight)*psikdxii*Fstar;
        residual = residual  - fterm;
        
        DFstardq = DFlux_dq_x(wstar,dpstardn);
        Jacobian = Jacobian - (nuv1*weight)*psikdxii*DFstardq*vectpsi(:,:,tquad,v1quad,v2quad,v3quad);

        %nortsout
        psikdeta = vectpsi_deta_Trans(:,:,tquad,v1quad,v2quad,v3quad);
        Gstar = Flux_y(wstar,pstar);
        gterm = (nuv2*weight)*psikdeta*Gstar;
        residual = residual  - gterm;
        
        DGstardq = DFlux_dq_y(wstar,dpstardn);
        Jacobian = Jacobian - (nuv2*weight)*psikdeta*DGstardq*vectpsi(:,:,tquad,v1quad,v2quad,v3quad);

        %upprdown
        psikdzta = vectpsi_dzta_Trans(:,:,tquad,v1quad,v2quad,v3quad);
        Hstar = Flux_z(wstar,pstar);
        gterm = (nuv3*weight)*psikdzta*Hstar;
        residual = residual  - gterm;
        
        DGstardq = DFlux_dq_z(wstar,dpstardn);
        Jacobian = Jacobian - (nuv3*weight)*psikdzta*DGstardq*vectpsi(:,:,tquad,v1quad,v2quad,v3quad);
        
       
        % source         
        sterm = data.deltat*Source(wstar)*weight;
        residual = residual - vectpsi_Trans(:,:,tquad,v1quad,v2quad,v3quad)*sterm;
        
        DstermdQ = data.deltat*DSourcedQ(wstar)*weight;
        Jacobian = Jacobian - vectpsi_Trans(:,:,tquad,v1quad,v2quad,v3quad)*DstermdQ*vectpsi(:,:,tquad,v1quad,v2quad,v3quad);
        
        
end
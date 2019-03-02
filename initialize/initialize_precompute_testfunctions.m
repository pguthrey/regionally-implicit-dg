function [data] = initialize_precompute_testfunctions(data)
% Precomputes the test functions at quadrature points
% written by Pierson Guthrey

Neqns = data.Neqns;
P = data.P;
Pt = P;
Px = P;

Py = 1;
Pz = 1;
if data.space_dims >= 2
    Py = P;
    if data.space_dims >= 3
        Pz = P;
    end
end

Pxplot = data.Pplot;
Pyplot = 1;
Pzplot = 1;
if data.space_dims >= 2
    Pyplot = data.Pplot;
    if data.space_dims >= 3
        Pzplot = data.Pplot;
    end
end


data.locs = data.locs;
thetaT = data.thetaT;
theta = data.theta;
sysBO = data.sysBO;
sysBOcoeffs = data.sysBOcoeffs;
sysBOs = data.sysBOs;
sysBOscoeffs = data.sysBOscoeffs;
 
vectpsi = zeros(Neqns,thetaT,Pt,Px,Py,Pz);
vectphi = zeros(Neqns,theta,Px,Py,Pz);
vectpsi_Trans = zeros(thetaT,Neqns,Pt,Px,Py,Pz);
vectphi_Trans = zeros(theta,Neqns,Px,Py,Pz);
vectphiplot = zeros(Neqns,theta,Pxplot,Pyplot,Pzplot);

syms xii real   
test = sym(zeros(1,data.M));
test_diff = sym(zeros(1,data.M));

testquad = zeros(1,data.M,length(data.locs));
for iell = 1:data.M
for iloc = 1:length(data.locs)
    testquad(1,iell,iloc) = testfunction_Legendre(data.locs(iloc),iell);    
end
end

plotquad = zeros(1,data.M,length(data.plotlocs));
for iell = 1:data.M
for iloc = 1:length(data.plotlocs)
    plotquad(1,iell,iloc) = testfunction_Legendre(data.plotlocs(iloc),iell);    
end
end


Zphi = 1;
Yphi = 1;

for kay = 1:theta
    for eqn = 1:Neqns
    ells = sysBOs(eqn,kay);
    coeffs = sysBOscoeffs(eqn,kay);
    for k = 1:data.Pd
        v1quad = data.Dlist(k,1); 
        if data.space_dims >= 2
            v2quad = data.Dlist(k,2);
        else
            v2quad = 1;
        end
        if data.space_dims >= 3
            v3quad = data.Dlist(k,3);
        else
            v3quad = 1;
        end
        if coeffs ~= 0 
            Xphi = testquad(1,data.BOs(ells,1),v1quad);
            if data.space_dims >= 2
                Yphi = testquad(1,data.BOs(ells,2),v2quad);
                if data.space_dims >= 3
                    Zphi = testquad(1,data.BOs(ells,3),v3quad);
                end
            end
        else
            Xphi = 0;
        end
        %Phi terms
        vectphi(eqn,kay,v1quad,v2quad,v3quad) =  Xphi.*Yphi.*Zphi.*coeffs;
        vectphi_Trans(kay,eqn,v1quad,v2quad,v3quad) =  vectphi(eqn,kay,v1quad,v2quad,v3quad);
    end
    end
end

for kay = 1:theta
    for eqn = 1:Neqns
    ells = sysBOs(eqn,kay);
    coeffs = sysBOscoeffs(eqn,kay);
    for k = 1:data.Pdplot
        v1quad = data.plotlist(k,1); 
        if data.space_dims >= 2
            v2quad = data.plotlist(k,2);
        else
            v2quad = 1;
        end
        if data.space_dims >= 3
            v3quad = data.plotlist(k,3);
        else
            v3quad = 1;
        end
        if coeffs ~= 0 
            Xphi = plotquad(1,data.BOs(ells,1),v1quad);
            if data.space_dims >= 2
                Yphi = plotquad(1,data.BOs(ells,2),v2quad);
                if data.space_dims >= 3
                    Zphi = plotquad(1,data.BOs(ells,3),v3quad);
                end
            end
        else
            Xphi = 0;
        end
        %Phi terms
        vectphiplot(eqn,kay,v1quad,v2quad,v3quad) =  Xphi.*Yphi.*Zphi.*coeffs;
    end
    end
end












Zpsi = 1;
Ypsi = 1;

for kay = 1:thetaT
        for eqn = 1:Neqns
        ell = sysBO(eqn,kay);
        coeffs = sysBOcoeffs(eqn,kay);
        for k = 1:data.Pdp1
            tquad = data.Dp1list(k,1); 
            v1quad = data.Dp1list(k,2); 
            if data.space_dims >= 2
                v2quad = data.Dp1list(k,3);
            else
                v2quad = 1;
            end
            if data.space_dims >= 3
                v3quad = data.Dp1list(k,4);
            else
                v3quad = 1;
            end 
        if coeffs ~= 0 
            Tpsi = testquad(1,data.BO(ell,1),tquad);
            Xpsi = testquad(1,data.BO(ell,2),v1quad);
            if data.space_dims >= 2
                Ypsi = testquad(1,data.BO(ell,3),v2quad);
                if data.space_dims >= 3
                    Zpsi = testquad(1,data.BO(ell,4),v3quad);
                end
            end
        end            
        vectpsi(eqn,kay,tquad,v1quad,v2quad,v3quad) =  Tpsi.*Xpsi.*Ypsi.*Zpsi.*coeffs;
        vectpsi_Trans(kay,eqn,tquad,v1quad,v2quad,v3quad) =  vectpsi(eqn,kay,tquad,v1quad,v2quad,v3quad);
        end
        end
end

alllocs = [-1 ; data.locs ; 1];
for kay = 1:thetaT
        for eqn = 1:Neqns
            ell = sysBO(eqn,kay);
            coeffs = sysBOcoeffs(eqn,kay);
            if coeffs ~=0
                for tquad = 1:Pt+2
                    for v1quad = 1:Px+2
                        v2quad = 1;
                        v3quad = 1;
                        tloc = alllocs(tquad);
                        v1loc = alllocs(v1quad);
                        Tpsi = testfunction_Legendre(alllocs(tquad),data.BO(ell,1));
                        Xpsi = testfunction_Legendre(alllocs(v1quad),data.BO(ell,2));       
                        vectpsi_all(eqn,kay,tquad,v1quad,v2quad,v3quad) =  Tpsi.*Xpsi.*coeffs;
                    end
                end
            end
        end
end

data.vectpsi_all = vectpsi_all;

data.vectpsi = vectpsi ;
data.vectpsi_Trans = vectpsi_Trans; 
data.vectphi = vectphi;
data.vectphi_Trans = vectphi_Trans;
data.vectphiplot = vectphiplot;

warning('set this value')
if ~strcmp(data.linearcase,'constantcoefficient')
  

    vectpsi_dxii = zeros(Neqns,thetaT,Pt,Px,Py,Pz);
    vectpsi_west = zeros(Neqns,thetaT,Px,Py,Pz);
    vectpsi_east = zeros(Neqns,thetaT,Px,Py,Pz);
    vectpsi_dxii_Trans = zeros(thetaT,Neqns,Px,Py,Pz);
    vectpsi_futr = zeros(Neqns,thetaT,Px,Py,Pz);
    vectpsi_dtau_Trans = zeros(thetaT,Neqns,Px,Py,Pz);
    vectpsi_futr_Trans = zeros(thetaT,Neqns,Px,Py,Pz);
    vectpsi_past_Trans = zeros(thetaT,Neqns,Px,Py,Pz);
    vectpsi_west_Trans = zeros(thetaT,Neqns,Px,Py,Pz);
    vectpsi_east_Trans = zeros(thetaT,Neqns,Px,Py,Pz);
    vectphi_dxii = zeros(Neqns,theta,Px,Py,Pz);
    vectphi_dxii_Trans = zeros(theta,Neqns,Px,Py,Pz);
    vectphi_west = zeros(Neqns,theta,Py,Pz);
    vectphi_east = zeros(Neqns,theta,Py,Pz);
    vectphi_west_Trans = zeros(theta,Neqns,Py,Pz);
    vectphi_east_Trans = zeros(theta,Neqns,Py,Pz);

    vectpsi_deta = zeros(Neqns,thetaT,Pt,Px,Py,Pz);
    vectpsi_nort = zeros(Neqns,thetaT,Px,Py,Pz);
    vectpsi_sout = zeros(Neqns,thetaT,Px,Py,Pz);
    vectpsi_deta_Trans = zeros(thetaT,Neqns,Px,Py,Pz);
    vectpsi_nort_Trans = zeros(thetaT,Neqns,Px,Py,Pz);
    vectpsi_sout_Trans = zeros(thetaT,Neqns,Px,Py,Pz);
    vectphi_deta = zeros(Neqns,theta,Px,Py,Pz);
    vectphi_deta_Trans = zeros(theta,Neqns,Px,Py,Pz);
    vectphi_nort = zeros(Neqns,theta,Py,Pz);
    vectphi_sout = zeros(Neqns,theta,Py,Pz);
    vectphi_nort_Trans = zeros(theta,Neqns,Py,Pz);
    vectphi_sout_Trans = zeros(theta,Neqns,Py,Pz);

    vectpsi_uppr = zeros(Neqns,thetaT,Px,Py,Pz);
    vectpsi_down = zeros(Neqns,thetaT,Px,Py,Pz);
    vectpsi_dzta = zeros(Neqns,thetaT,Pt,Px,Py,Pz);
    vectpsi_dzta_Trans = zeros(thetaT,Neqns,Px,Py,Pz);
    vectpsi_uppr_Trans = zeros(thetaT,Neqns,Px,Py,Pz);
    vectpsi_down_Trans = zeros(thetaT,Neqns,Px,Py,Pz);
    vectphi_dzta = zeros(Neqns,theta,Px,Py,Pz);
    vectphi_dzta_Trans = zeros(theta,Neqns,Px,Py,Pz);
    vectphi_uppr = zeros(Neqns,theta,Py,Pz);
    vectphi_down = zeros(Neqns,theta,Py,Pz);
    vectphi_uppr_Trans = zeros(theta,Neqns,Py,Pz);
    vectphi_down_Trans = zeros(theta,Neqns,Py,Pz);

    for kay = 1:thetaT
        for eqn = 1:Neqns
        ell = sysBO(eqn,kay);
        coeff = sysBOcoeffs(eqn,kay);

        for k = 1:data.Pd
            v1quad = data.Dlist(k,1); 
            v1loc = data.Dquadlocs(k,1);
            if data.space_dims >= 2
                v2quad = data.Dlist(k,2);
                v2loc = data.Dquadlocs(k,2);
            else
                v2quad = 1;
                v2loc = inf;
            end
            if data.space_dims >= 3
                v3quad = data.Dlist(k,3);
                v3loc = data.Dquadlocs(k,3);
            else
                v3quad = 1;
                v3loc = inf;
            end
            vectpsi_west(eqn,kay,v1quad,v2quad,v3quad) = testfunction_psi([v1loc,-1,v2loc,v3loc],ell,data).*coeff;
            vectpsi_east(eqn,kay,v1quad,v2quad,v3quad) = testfunction_psi([v1loc, 1,v2loc,v3loc],ell,data).*coeff;
            vectpsi_west_Trans(kay,eqn,v1quad,v2quad,v3quad) = testfunction_psi([v1loc,-1,v2loc,v3loc],ell,data).*coeff;
            vectpsi_east_Trans(kay,eqn,v1quad,v2quad,v3quad) = testfunction_psi([v1loc, 1,v2loc,v3loc],ell,data).*coeff;       
            if data.space_dims >= 2
                vectpsi_nort(eqn,kay,v1quad,v2quad,v3quad) = testfunction_psi([v1loc,v2loc, 1,v3loc],ell,data).*coeff;
                vectpsi_sout(eqn,kay,v1quad,v2quad,v3quad) = testfunction_psi([v1loc,v2loc,-1,v3loc],ell,data).*coeff;
                vectpsi_nort_Trans(kay,eqn,v1quad,v2quad,v3quad) = testfunction_psi([v1loc,v2loc, 1,v3loc],ell,data).*coeff;
                vectpsi_sout_Trans(kay,eqn,v1quad,v2quad,v3quad) = testfunction_psi([v1loc,v2loc,-1,v3loc],ell,data).*coeff;            
                if data.space_dims >= 3
                    vectpsi_uppr(eqn,kay,v1quad,v2quad,v3quad) = testfunction_psi([v1loc,v2loc,v3loc,1],ell,data).*coeff;
                    vectpsi_down(eqn,kay,v1quad,v2quad,v3quad) = testfunction_psi([v1loc,v2loc,v3loc,-1],ell,data).*coeff;
                    vectpsi_uppr_Trans(kay,eqn,v1quad,v2quad,v3quad) = testfunction_psi([v1loc,v2loc,v3loc,1],ell,data).*coeff;
                    vectpsi_down_Trans(kay,eqn,v1quad,v2quad,v3quad) = testfunction_psi([v1loc,v2loc,v3loc,-1],ell,data).*coeff;
                end
            end
            vectpsi_futr(eqn,kay,v1quad,v2quad,v3quad) =  testfunction_psi([1,v1loc,v2loc,v3loc],ell,data).*coeff;
            vectpsi_futr_Trans(kay,eqn,v1quad,v2quad,v3quad) =  testfunction_psi([1,v1loc,v2loc,v3loc],ell,data).*coeff;
            vectpsi_past_Trans(kay,eqn,v1quad,v2quad,v3quad) =  testfunction_psi([-1,v1loc,v2loc,v3loc],ell,data).*coeff;
        end
        end
    end

    clear v1quad v2quad v3quad v1loc v2loc v3loc tquad tloc ell coeff ells coeffs quadpoint

    %Phi Terms
    for kay = 1:theta
        for eqn = 1:Neqns
        ells = sysBOs(eqn,kay);
        coeffs = sysBOscoeffs(eqn,kay);
        for k = 1:data.Pdm1
            if data.space_dims >= 2
                v1quad = data.Dm1list(k,1);
                v1loc = data.Dm1quadlocs(k,1);
            else
                v1quad = 1;
                v1loc = inf;
            end
            if data.space_dims >= 3
                v2quad = data.Dm1list(k,2);
                v2loc = data.Dm1quadlocs(k,2);
            else
                v2quad = 1;
                v2loc = inf;
            end
            vectphi_east(eqn,kay,v1quad,v2quad) =  testfunction_phi([1,v1loc,v2loc],ells,data).*coeffs;
            vectphi_east_Trans(kay,eqn,v1quad,v2quad) =  testfunction_phi([1,v1loc,v2loc],ells,data).*coeffs;        
            vectphi_west(eqn,kay,v1quad,v2quad) =  testfunction_phi([-1,v1loc,v2loc],ells,data).*coeffs;
            vectphi_west_Trans(kay,eqn,v1quad,v2quad) =  testfunction_phi([-1,v1loc,v2loc],ells,data).*coeffs;
            if data.space_dims >= 2
                vectphi_nort(eqn,kay,v1quad,v2quad) =  testfunction_phi([v1loc,1,v2loc],ells,data).*coeffs;
                vectphi_nort_Trans(kay,eqn,v1quad,v2quad) =  testfunction_phi([v1loc,1,v2loc],ells,data).*coeffs;        
                vectphi_sout(eqn,kay,v1quad,v2quad) =  testfunction_phi([v1loc,-1,v2loc],ells,data).*coeffs;
                vectphi_sout_Trans(kay,eqn,v1quad,v2quad) =  testfunction_phi([v1loc,-1,v2loc],ells,data).*coeffs;
                if data.space_dims >= 3
                    vectphi_uppr(eqn,kay,v1quad,v2quad) =  testfunction_phi([v1loc,v2loc,1],ells,data).*coeffs;
                    vectphi_uppr_Trans(kay,eqn,v1quad,v2quad) =  testfunction_phi([v1loc,v2loc,1],ells,data).*coeffs;        
                    vectphi_down(eqn,kay,v1quad,v2quad) =  testfunction_phi([v1loc,v2loc,-1],ells,data).*coeffs;
                    vectphi_down_Trans(kay,eqn,v1quad,v2quad) =  testfunction_phi([v1loc,v2loc,-1],ells,data).*coeffs;
                end
            end
        end
        end
    end

    clear v1quad v2quad v3quad v1loc v2loc v3loc tquad tloc ell coeff ells coeffs quadpoint

    for kay = 1:theta
        for eqn = 1:Neqns
        ells = sysBOs(eqn,kay);
        coeffs = sysBOscoeffs(eqn,kay);
        for k = 1:data.Pd
            v1quad = data.Dlist(k,1); 
            v1loc = data.Dquadlocs(k,1);
            if data.space_dims >= 2
                v2quad = data.Dlist(k,2);
                v2loc = data.Dquadlocs(k,2);
            else
                v2quad = 1;
                v2loc = inf;
            end
            if data.space_dims >= 3
                v3quad = data.Dlist(k,3);
                v3loc = data.Dquadlocs(k,3);
            else
                v3quad = 1;
                v3loc = inf;
            end
            quadpoint = [v1loc,v2loc,v3loc];
            %Phi terms
            vectphi(eqn,kay,v1quad,v2quad,v3quad) =  testfunction_phi(quadpoint,ells,data).*coeffs;
            vectphi_Trans(kay,eqn,v1quad,v2quad,v3quad) =  testfunction_phi(quadpoint,ells,data).*coeffs;
            vectphi_dxii(eqn,kay,v1quad,v2quad,v3quad) = testfunction_phi_dxii(quadpoint,ells,data).*coeffs;
            if data.space_dims >= 2
                vectphi_deta(eqn,kay,v1quad,v2quad,v3quad) = testfunction_phi_deta(quadpoint,ells,data).*coeffs;
                if data.space_dims >= 3
                    vectphi_dzta(eqn,kay,v1quad,v2quad,v3quad) = testfunction_phi_dzta(quadpoint,ells,data).*coeffs;
                end
            end
        end
        end
    end
 
    clear v1quad v2quad v3quad v1loc v2loc v3loc tquad tloc ell coeff ells coeffs quadpoint


    %Phi derivative terms
    for kay = 1:theta
        for eqn = 1:Neqns
        ells = sysBOs(eqn,kay);
        coeffs = sysBOscoeffs(eqn,kay);
        for k = 1:data.Pd
            v1quad = data.Dlist(k,1); 
            v1loc = data.Dquadlocs(k,1);
            if data.space_dims >= 2
                v2quad = data.Dlist(k,2);
                v2loc = data.Dquadlocs(k,2);
            else
                v2quad = 1;
                v2loc = inf;
            end
            if data.space_dims >= 3
                v3quad = data.Dlist(k,3);
                v3loc = data.Dquadlocs(k,3);
            else
                v3quad = 1;
                v3loc = inf;
            end
            quadpoint = [v1loc,v2loc,v3loc];
            %Phi terms
            vectphi_dxii(eqn,kay,v1quad,v2quad,v3quad) = testfunction_phi_dxii(quadpoint,ells,data).*coeffs;
            vectphi_dxii_Trans(kay,eqn,v1quad,v2quad,v3quad) = testfunction_phi_dxii(quadpoint,ells,data).*coeffs;
            if data.space_dims >= 2
                vectphi_deta(eqn,kay,v1quad,v2quad,v3quad) = testfunction_phi_deta(quadpoint,ells,data).*coeffs;
                vectphi_deta_Trans(kay,eqn,v1quad,v2quad,v3quad) = testfunction_phi_deta(quadpoint,ells,data).*coeffs;
                if data.space_dims >= 3
                    vectphi_dzta(eqn,kay,v1quad,v2quad,v3quad) = testfunction_phi_dzta(quadpoint,ells,data).*coeffs;
                    vectphi_dzta_Trans(kay,eqn,v1quad,v2quad,v3quad) = testfunction_phi_dzta(quadpoint,ells,data).*coeffs;
                end
            end
        end
        end
    end

    clear v1quad v2quad v3quad v1loc v2loc v3loc tquad tloc ell coeff ells coeffs quadpoint

    for kay = 1:thetaT
        for eqn = 1:Neqns
        ell = sysBO(eqn,kay);
        coeff = sysBOcoeffs(eqn,kay);
        for k = 1:data.Pdp1
            tquad = data.Dp1list(k,1); 
            tloc = data.Dp1quadlocs(k,1);
            v1quad = data.Dp1list(k,2); 
            v1loc = data.Dp1quadlocs(k,2);
            if data.space_dims >= 2
                v2quad = data.Dp1list(k,3);
                v2loc = data.Dp1quadlocs(k,3);
            else
                v2quad = 1;
                v2loc = inf;
            end
            if data.space_dims >= 3
                v3quad = data.Dp1list(k,4);
                v3loc = data.Dp1quadlocs(k,4);
            else
                v3quad = 1;
                v3loc = inf;
            end
            quadpoint = [tloc,v1loc,v2loc,v3loc];
            vectpsi_dtau_Trans(kay,eqn,tquad,v1quad,v2quad,v3quad) = testfunction_psi_dtau(quadpoint,ell,data).*coeff;
            vectpsi_dxii(eqn,kay,tquad,v1quad,v2quad,v3quad) = testfunction_psi_dxii(quadpoint,ell,data).*coeff;
            vectpsi_dxii_Trans(kay,eqn,tquad,v1quad,v2quad,v3quad) = testfunction_psi_dxii(quadpoint,ell,data).*coeff;
            if data.space_dims >= 2
                vectpsi_deta_Trans(kay,eqn,tquad,v1quad,v2quad,v3quad) = testfunction_psi_deta(quadpoint,ell,data).*coeff;         
                vectpsi_deta(eqn,kay,tquad,v1quad,v2quad,v3quad) = testfunction_psi_deta(quadpoint,ell,data).*coeff;
                if data.space_dims >= 3
                vectpsi_dzta(eqn,kay,tquad,v1quad,v2quad,v3quad) = testfunction_psi_dzta(quadpoint,ell,data).*coeff;
                vectpsi_dzta_Trans(kay,eqn,tquad,v1quad,v2quad,v3quad) = testfunction_psi_dzta(quadpoint,ell,data).*coeff;
                end
            end        
        end
        end
    end


    clear v1quad v2quad v3quad v1loc v2loc v3loc tquad tloc ell coeff ells coeffs quadpoint

    for kay = 1:theta
        for eqn = 1:Neqns
        ells = sysBOs(eqn,kay);
        coeffs = sysBOscoeffs(eqn,kay);
        for k = 1:data.Pdlimiter
            v1quad = data.limiterlist(k,1); 
            v1loc = data.Dlimiterlocs(k,1);
            if data.space_dims >= 2
                v2quad = data.limiterlist(k,2);
                v2loc = data.Dlimiterlocs(k,2);
            else
                v2quad = 1;
                v2loc = inf;
            end
            if data.space_dims >= 3
                v3quad = data.limiterlist(k,3);
                v3loc = data.Dlimiterlocs(k,3);
            else
                v3quad = 1;
                v3loc = inf;
            end
            quadpoint = [v1loc,v2loc,v3loc];
            %Phi terms
            data.vectphi_limiter(eqn,kay,v1quad,v2quad,v3quad) =  testfunction_phi(quadpoint,ells,data).*coeffs;        
        end
        end
    end

    for kay = 1:thetaT
        for eqn = 1:Neqns
        ells = sysBO(eqn,kay);
        coeffs = sysBOcoeffs(eqn,kay);
        for tquad = 1:data.P+2
        for v1quad = 1:data.P+2
        for v2quad = 1:data.P+2        
            tloc = data.limiterlocs(tquad);
            v1loc = data.limiterlocs(v1quad);
            v2loc = data.limiterlocs(v2quad);
            v3loc = inf;
            quadpoint = [tloc,v1loc,v2loc,v3loc];
            data.vectpsi_limiter(eqn,kay,tquad,v1quad,v2quad,v3quad) =  testfunction_psi(quadpoint,ells,data).*coeffs;        
        end
        end
        end
        end
    end
    %{
    limiter_locs = [-1;data.locs;1];
    data.vectphi_limiter = zeros(Neqns,theta,P+2,P+2);
    for kay = 1:theta
    for eqn = 1:Neqns
        ells = sysBOs(eqn,kay);
        coeffs = sysBOscoeffs(eqn,kay);   
        for v1quad = 1:(P+2)
        for v2quad = 1:(P+2)
        for v3quad = 1:(P+2)
            xii = limiter_locs(v1quad);
            eta = limiter_locs(v2quad);
            zta = limiter_locs(v3quad);
            data.vectphi_limiter(eqn,kay,v1quad,v2quad,v3quad) =  testfunction_phi([xii,eta,zta],ells,data).*coeffs;
        end
        end
        end
    end
    end
    %}
    %{
    data.vectpsi_limiter = zeros(Neqns,thetaT,P+2,P+2);
    for kay = 1:thetaT
    for eqn = 1:Neqns
        ells = sysBO(eqn,kay);
        coeffs = sysBOcoeffs(eqn,kay);   
        for v1quad = 1:(P+2)
        for v2quad = 1:(P+2)
        for v3quad = 1:(P+2)
        for tquad = 1:(P+2)
            xii = limiter_locs(v1quad);
            eta = limiter_locs(v2quad);
            zta = limiter_locs(v3quad);
            tau = limiter_locs(tquad);
            data.vectphi_limiter(eqn,kay,v1quad,v2quad,v3quad) =  testfunction_phi([xii,eta,zta],ells,data).*coeffs;
        end
        end
        end
        end
    end
    end
    %}
    %}

    %%{


    data.vectphi_average = ((data.sysBOs == 1).*data.sysBOscoeffs);
    data.vectpsi_average = ((data.sysBO == 1).*data.sysBOcoeffs);
    if ~strcmp(data.basis,'canonical')
        warning('Phi and psi averages possibly not computed properly')
        warning('int of  phi^T phi possibly not computer properly')
    end

    data.vectpsi_dxii = vectpsi_dxii ;
    data.vectpsi_futr = vectpsi_futr;
    data.vectpsi_west = vectpsi_west;
    data.vectpsi_east = vectpsi_east;
    data.vectpsi_dxii_Trans = vectpsi_dxii_Trans;
    data.vectpsi_futr_Trans = vectpsi_futr_Trans;
    data.vectpsi_past_Trans = vectpsi_past_Trans;
    data.vectpsi_west_Trans = vectpsi_west_Trans;
    data.vectpsi_east_Trans = vectpsi_east_Trans;
    data.vectphi_dxii = vectphi_dxii;
    data.vectphi_dxii_Trans = vectphi_dxii_Trans;
    data.vectphi_west = vectphi_west;
    data.vectphi_east = vectphi_east;
    data.vectphi_west_Trans = vectphi_west_Trans;
    data.vectphi_east_Trans = vectphi_east_Trans;
    data.vectpsi_dtau_Trans = vectpsi_dtau_Trans;

    if data.space_dims >= 2 
        data.vectpsi_deta = vectpsi_deta;
        data.vectpsi_nort = vectpsi_nort;
        data.vectpsi_sout = vectpsi_sout;
        data.vectpsi_deta_Trans = vectpsi_deta_Trans;
        data.vectpsi_nort_Trans = vectpsi_nort_Trans;
        data.vectpsi_sout_Trans = vectpsi_sout_Trans;
        data.vectphi_nort = vectphi_nort;
        data.vectphi_sout = vectphi_sout;
        data.vectphi_nort_Trans = vectphi_nort_Trans;
        data.vectphi_sout_Trans = vectphi_sout_Trans;
        data.vectphi_deta = vectphi_deta;    
        data.vectphi_deta_Trans = vectphi_deta_Trans;
        if data.space_dims >= 3
            data.vectpsi_dzta = vectpsi_dzta;
            data.vectpsi_uppr = vectpsi_uppr;
            data.vectpsi_down = vectpsi_down;
            data.vectpsi_dzta_Trans = vectpsi_dzta_Trans;
            data.vectpsi_uppr_Trans = vectpsi_uppr_Trans;
            data.vectpsi_down_Trans = vectpsi_down_Trans;
            data.vectphi_dzta = vectphi_dzta;
            data.vectphi_dzta_Trans = vectphi_dzta_Trans;
            data.vectphi_uppr = vectphi_uppr;
            data.vectphi_down = vectphi_down;
            data.vectphi_uppr_Trans = vectphi_uppr_Trans;
            data.vectphi_down_Trans = vectphi_down_Trans;
        end
    end

    data.bases = zeros(data.thetaT,data.theta,Pt,Px,Py,Pz);

    if strcmp(data.methodtype,'explicit')
        locs = data.locs;
        torder = data.nstages/2+1;
        Bmat = data.predictorB;
        %Stuff for the explicit method 
        for k = 1:data.Pdp1
            tquad = data.Dp1list(k,1);
            v1quad = data.Dp1list(k,2); 
            if data.space_dims >= 2
                v2quad = data.Dp1list(k,3);
            else
                v2quad = 1;
            end
            if data.space_dims >= 3
                v3quad = data.Dp1list(k,4);
            else
                v3quad = 1;
            end
            psiT = data.vectpsi_Trans(:,:,tquad,v1quad,v2quad,v3quad);
            phi = data.vectphi(:,:,v1quad,v2quad,v3quad);       
            data.bases(:,:,tquad,v1quad,v2quad,v3quad) = psiT*phi;

            tloc = (locs(tquad)+1)/2;
            Tvect = ((tloc).^((torder:-1:1)))';        
            data.Cmat(:,:,tquad,v1quad,v2quad,v3quad) = Bmat*Tvect;

        end  
    end


    for kay = 1:thetaT
        for eqn = 1:Neqns
        ells = sysBO(eqn,kay);
        coeffs = sysBOcoeffs(eqn,kay);
        for k = 1:data.Pdplotpsi
            tquad = data.plotlistpsi(k,1); 
            tloc = data.plotlistpsi(k,1); 
            v1quad = data.plotlistpsi(k,2); 
            v1loc = data.Dplotlocspsi(k,2);
            if data.space_dims >= 2
                v2quad = data.plotlistpsi(k,3);
                v2loc = data.Dplotlocspsi(k,3);
            else
                v2quad = 1;
                v2loc = inf;
            end
            if data.space_dims >= 3
                v3quad = data.plotlistpsi(k,4);
                v3loc = data.Dplotlocspsi(k,4);
            else
                v3quad = 1;
                v3loc = inf;
            end
            quadpoint = [tloc,v1loc,v2loc,v3loc];
            %Phi terms
            data.vectpsiplot(eqn,kay,tquad,v1quad,v2quad,v3quad) =  testfunction_psi(quadpoint,ells,data).*coeffs;        
        end
        end
    end

    data.inner_prod_average_phi = data.vectphi_average'*data.vectphi_average;
    data.inner_prod_average_psi = data.vectpsi_average'*data.vectpsi_average;

    if data.space_dims == 2

    for kay = 1:theta
        for eqn = 1:Neqns
        ells = sysBOs(eqn,kay);
        coeffs = sysBOscoeffs(eqn,kay);
        for k = 1:length(data.plotlocs)
            v1loc = data.plotlocs(k);
            vectphi_2Dto1D(eqn,kay,k) =  testfunction_phi([v1loc,0,v2loc],ells,data).*coeffs;
        end
        end
    end
        data.vectphi_2Dto1D = vectphi_2Dto1D;    
    end
end

v1quad = 1;
v2quad = 1;
v3quad = 1;

for k = 1:data.Pd
    v1quad = data.Dlist(k,1); 
    if data.space_dims >= 2
        v2quad = data.Dlist(k,2);
        if data.space_dims >= 3
            v3quad = data.Dlist(k,3);
        end        
    end

    psie = data.vectpsi_east(:,:,v1quad,v2quad,v3quad); 
    psiw = data.vectpsi_west(:,:,v1quad,v2quad,v3quad); 

    data.psikpsim_ee(:,:,v1quad,v2quad,v3quad) = psie'*psie; 
    data.psikpsim_we(:,:,v1quad,v2quad,v3quad) = psiw'*psie; 
    data.psikpsim_ew(:,:,v1quad,v2quad,v3quad) = psie'*psiw; 
    data.psikpsim_ww(:,:,v1quad,v2quad,v3quad) = psiw'*psiw; 
    
    if data.space_dims >= 2
        psin = data.vectpsi_nort(:,:,v1quad,v2quad,v3quad); 
        psis = data.vectpsi_sout(:,:,v1quad,v2quad,v3quad); 

        data.psikpsim_nn(:,:,v1quad,v2quad,v3quad) = psin'*psin; 
        data.psikpsim_sn(:,:,v1quad,v2quad,v3quad) = psis'*psin; 
        data.psikpsim_ns(:,:,v1quad,v2quad,v3quad) = psin'*psis; 
        data.psikpsim_ss(:,:,v1quad,v2quad,v3quad) = psis'*psis;     
    end
end
    





%{
syms t real
syms x real
vectpsi = sym(zeros(thetaT,Neqns));
vectphi = sym(zeros(theta,Neqns));
for kay = 1:thetaT
    for eqn = 1:Neqns
        ell = sysBO(eqn,kay);
        coeffs = sysBOcoeffs(eqn,kay);
        vectpsi(kay,eqn) = coeffs.*testfunction_psi([t,x,NaN,NaN],ell,data);
    end
end
for kay = 1:theta
    for eqn = 1:Neqns
        ell = sysBOs(eqn,kay);
        coeffs = sysBOscoeffs(eqn,kay);
        vectphi(kay,eqn) = testfunction_phi([x,NaN,NaN],ell,data);
    end
end
Lmat = int(int(vectpsi*diff(vectpsi',t),t,-1,1),x,-1,1)/4 ...
    + int(subs(vectpsi*vectpsi',t,-1),x,-1,1)/4;
data.Lmat = double(Lmat);
backmat = Lmat\int(subs(vectpsi,t,-1)*vectphi',x,-1,1)/4;
data.backmat = double(backmat);
data.Tnowinv = inv(data.Lmat)/4;
%}








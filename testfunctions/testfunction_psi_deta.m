function [ psi ] = testfunction_psi_deta(quadpoint,ell,data)
% Evalautes partial_eta psi^(ell)(tau,xii,eta,zeta)
% written by Pierson Guthrey


psi_zta = 1;
psi_eta = 1;

if ell ~=0
    ell_tau = data.BO(ell,1);
    tau = quadpoint(1);
    psi_tau = testfunction_Legendre(tau,ell_tau);
    ell_xii = data.BO(ell,2);
    xii = quadpoint(2);
    psi_xii = testfunction_Legendre(xii,ell_xii);
    if data.space_dims >= 2
        eta = quadpoint(3);
        ell_eta = data.BO(ell,3);
        psi_eta = testfunction_LegendreDiff(eta,ell_eta);
        if data.space_dims >= 3
            zta = quadpoint(4);
            ell_zta = data.BO(ell,4);
            psi_zta = testfunction_Legendre(zta,ell_zta);
        end
    end           
    psi = psi_tau.*psi_xii.*psi_eta.*psi_zta;
else
    psi = 0;    
end
end
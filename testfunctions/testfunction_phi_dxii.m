function [ phi ] = testfunction_phi_dxii(quadpoint,ell,data)
% Evalautes partail_xii phi^(ell)(xii,eta,zeta)
% written by Pierson Guthrey

phi_zta = 1;
phi_eta = 1;

if ell ~=0
    ell_xii = data.BOs(ell,1);
    xii = quadpoint(1);
    phi_xii = testfunction_LegendreDiff(xii,ell_xii);
    if data.space_dims >= 2
        eta = quadpoint(2);
        ell_eta = data.BOs(ell,2);
        phi_eta = testfunction_Legendre(eta,ell_eta);
        if data.space_dims >= 3
            zta = quadpoint(3);
            ell_zta = data.BOs(ell,3);
            phi_zta = testfunction_Legendre(zta,ell_zta);
        end
    end           
    phi = phi_xii.*phi_eta.*phi_zta;
else
    phi = 0;
end
end
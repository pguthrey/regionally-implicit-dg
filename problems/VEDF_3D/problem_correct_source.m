function [correct_cell,auxiliary,data] = problem_correct_source(correct_cell,qstar,auxiliary,data,cellcenter,cellindex)

vectphi_Trans = data.vectphi_Trans;
vectpsi = data.vectpsi;

%{
Pdp1 = data.Pdp1;
space_dims = data.space_dims;
Dp1list = data.Dp1list;
Dp1quadwgts = data.Dp1quadwgts;
%}

vectphi = data.vectphi;
Pdp1 = data.Pdp1;
space_dims = data.space_dims;
Dp1list = data.Dp1list;
Dp1quadwgts = data.Dp1quadwgts;

tau = 1;

v2quad = 1;
v3quad = 1;

total_viscocity = auxiliary.total_viscocity(:,cellindex(1),cellindex(2),cellindex(3));

% EXTERNAL FORCES SET TO ZERO
Source = @(q)[  0;
                q(5)/q(1);
                q(6)/q(1);
                q(7)/q(1);
                (total_viscocity(1)-q(5))/tau
                (total_viscocity(2)-q(6))/tau
                (total_viscocity(3)-q(7))/tau ];

for k = 1:Pdp1
        tquad = Dp1list(k,1);
        v1quad = Dp1list(k,2);
        if space_dims >= 2
            v2quad = Dp1list(k,3);
            if space_dims >= 3
                v3quad = Dp1list(k,4);
            end
        end
        weight = Dp1quadwgts(k);
        wstar = vectpsi(:,:,tquad,v1quad,v2quad,v3quad)*qstar;
        
        source = data.deltat*Source(wstar);
        correct_cell = correct_cell - vectphi_Trans(:,:,v1quad,v2quad,v3quad)*source*weight;
end






end


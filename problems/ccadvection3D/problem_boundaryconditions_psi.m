function [DGprev_east,DGprev_west,DGprev_nort,DGprev_sout] = problem_boundaryconditions_psi(DGprev,data)


error('whoops')
DGprev_east = NaN(data.thetaT,data.Nv1,data.Nv2,data.Nv3);
DGprev_west = NaN(data.thetaT,data.Nv1,data.Nv2,data.Nv3);
DGprev_nort = NaN(data.thetaT,data.Nv1,data.Nv2,data.Nv3);
DGprev_sout = NaN(data.thetaT,data.Nv1,data.Nv2,data.Nv3);
DGprev_uppr = NaN(data.thetaT,data.Nv1,data.Nv2,data.Nv3);
DGprev_down = NaN(data.thetaT,data.Nv1,data.Nv2,data.Nv3);

for iv2 = 1:data.Nv2
for iv3 = 1:data.Nv3
    DGprev_east(:,:,iv2,iv3) = DGprev(:,[2:end 1],iv2,iv3);
    DGprev_west(:,:,iv2,iv3) = DGprev(:,[end 1:end-1],iv2,iv3);
end
end

for iv1 = 1:data.Nv1
for iv3 = 1:data.Nv3
    DGprev_nort(:,iv1,:,iv3) = DGprev(:,iv1,[2:end 1],iv3);
    DGprev_sout(:,iv1,:,iv3) = DGprev(:,iv1,[end 1:end-1],iv3);
end
end

for iv1 = 1:data.Nv1
for iv2 = 1:data.Nv2
    DGprev_uppr(:,iv1,iv2,:) = DGprev(:,iv1,iv2,[2:end 1]);
    DGprev_down(:,iv1,iv2,:) = DGprev(:,iv1,iv2,[end 1:end-1]);
end
end



end


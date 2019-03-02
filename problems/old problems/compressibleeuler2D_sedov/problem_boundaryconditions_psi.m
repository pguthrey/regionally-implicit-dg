function [DGprev_east,DGprev_west,DGprev_nort,DGprev_sout] = problem_boundaryconditions_psi(DGprev,data)

DGprev_east = NaN(data.thetaT,data.Nv1,data.Nv2);
DGprev_west = NaN(data.thetaT,data.Nv1,data.Nv2);
DGprev_nort = NaN(data.thetaT,data.Nv1,data.Nv2);
DGprev_sout = NaN(data.thetaT,data.Nv1,data.Nv2);


DGprev_east(:,1:(end-1),:) = DGprev(:,2:end,:);
for iv2= 1:data.Nv2
    DGprev_east(:,end,iv2) = data.extrapolate_psi_x*DGprev(:,end,iv2);
end

DGprev_west(:,2:(end),:) = DGprev(:,1:(end-1),:);
for iv2= 1:data.Nv2
    DGprev_west(:,1,iv2) = data.extrapolate_psi_x*DGprev(:,1,iv2);
end

DGprev_nort(:,:,1:(end-1)) = DGprev(:,:,2:end);
for iv1= 1:data.Nv1
    DGprev_nort(:,iv1,end) = data.extrapolate_psi_y*DGprev(:,iv1,end);
end

DGprev_sout(:,:,2:(end)) = DGprev(:,:,1:(end-1));
for iv1= 1:data.Nv1
    DGprev_sout(:,iv1,1) = data.extrapolate_psi_y*DGprev(:,iv1,1);
end


end


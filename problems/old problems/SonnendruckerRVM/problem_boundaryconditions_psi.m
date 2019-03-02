function [DGprev_east,DGprev_west,DGprev_nort,DGprev_sout] = problem_boundaryconditions_psi(DGprev,data)

DGprev_east = NaN(data.thetaT,data.Nv1,data.Nv2);
DGprev_west = NaN(data.thetaT,data.Nv1,data.Nv2);
DGprev_nort = NaN(data.thetaT,data.Nv1,data.Nv2);
DGprev_sout = NaN(data.thetaT,data.Nv1,data.Nv2);

DGprev_east = [DGprev(:,[2:end],:) 0.*DGprev(:,1,:)];
DGprev_west = [0.*DGprev(:,end,:) DGprev(:,1:(end-1),:)];

for iv1 = 1:data.Nv1
    DGprev_sout(:,iv1,2:end) = DGprev(:,iv1,1:(end-1));
    DGprev_sout(:,iv1,1) = 0.*DGprev(:,iv1,end);

    DGprev_nort(:,iv1,1:(end-1)) = DGprev(:,iv1,2:end);
    DGprev_nort(:,iv1,end) = 0.*DGprev(:,iv1,1);
end

end


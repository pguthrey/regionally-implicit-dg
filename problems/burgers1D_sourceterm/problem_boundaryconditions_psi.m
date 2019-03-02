function [DGprev_east,DGprev_west] = problem_boundaryconditions_psi(DGprev,data)

%{
DGprev_east = [DGprev(:,[2:end],:) data.extrapolate_psi_x*DGprev(:,end,:)];
DGprev_west = [data.extrapolate_psi_x*DGprev(:,1,:) DGprev(:,1:(end-1),:)];
%}
DGprev_east = [DGprev(:,[2:end],:) DGprev(:,1,:)];
DGprev_west = [DGprev(:,end,:) DGprev(:,1:(end-1),:)];

end


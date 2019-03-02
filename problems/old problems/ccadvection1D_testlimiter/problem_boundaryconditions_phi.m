function [DGprev_east,DGprev_west] = problem_boundaryconditions_phi(DGprev,data)

DGprev_east = [DGprev(:,2:end) 0.*DGprev(:,1)];
DGprev_west = [0.*DGprev(:,1) DGprev(:,1:(end-1))];
end


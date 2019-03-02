function [DGprev_east,DGprev_west] = problem_boundaryconditions_psi(DGprev,data)

DGprev_east = DGprev(:,[2:end 1]);
DGprev_west = DGprev(:,[end 1:(end-1)]);

end


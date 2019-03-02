function [DGprev_east,DGprev_west] = problem_boundaryconditions_phi(DGprev,data)

error('WHOOPS')

DGprev_east = DGprev(:,[2:end end]);
DGprev_west = DGprev(:,[1 1:(end-1)]);

end


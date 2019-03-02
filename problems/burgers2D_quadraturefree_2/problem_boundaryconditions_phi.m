function [DGprev_east,DGprev_west,DGprev_nort,DGprev_sout] = problem_boundaryconditions_phi(DGprev,data)

[DGprev_east,DGprev_west,DGprev_nort,DGprev_sout] = problem_boundaryconditions_psi(DGprev,data);

end


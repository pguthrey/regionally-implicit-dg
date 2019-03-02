function [DGprev_east,DGprev_west,DGprev_nort,DGprev_sout] = problem_boundaryconditions_phi(DGprev,data)

DGprev_east = [DGprev(:,[2:end],:) extrapolate_phi_x*DGprev(:,end,:)];
DGprev_west = [extrapolate_phi_x*DGprev(:,1,:) DGprev(:,1:(end-1),:)];

DGprev_nort = [extrapolate_phi_y*DGprev(:,:,1) ; DGprev(:,:,1:(end-1))];
DGprev_sout = [DGprev(:,:,2:end); extrapolate_phi_y*DGprev(:,:,end)];

end


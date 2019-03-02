function [DGprev_east,DGprev_west,DGprev_nort,DGprev_sout] = problem_boundaryconditions_psi(DGprev,data)


DGprev_east = 0.*DGprev;
DGprev_west = 0.*DGprev;
DGprev_nort = 0.*DGprev;
DGprev_sout = 0.*DGprev;

DGprev_east(:,1:(end-1),:) = DGprev(:,2:end,:) ;
DGprev_east(:,end,:) = DGprev(:,1,:);

DGprev_west(:,2:end,:) = DGprev(:,1:(end-1),:) ;
DGprev_west(:,1,:) = DGprev(:,end,:);

DGprev_nort(:,:,1:(end-1)) = DGprev(:,:,2:end) ;
DGprev_nort(:,:,end) = DGprev(:,:,1);

DGprev_sout(:,:,2:end) = DGprev(:,:,1:(end-1)) ;
DGprev_sout(:,:,1) = DGprev(:,:,end);

end


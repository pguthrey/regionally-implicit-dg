function [DGprev_east,auxiliary,DGprev_west,DGprev_nort,DGprev_sout,DGprev_uppr,DGprev_down] = problem_boundaryconditions_psi(DGprev,auxiliary,data)


DGprev_east = 0.*DGprev;
DGprev_west = 0.*DGprev;
DGprev_nort = 0.*DGprev;
DGprev_sout = 0.*DGprev;
DGprev_uppr = 0.*DGprev;
DGprev_down = 0.*DGprev;

DGprev_east(:,1:(end-1),:,:) = DGprev(:,2:end,:,:) ;
DGprev_east(:,end,:,:) = DGprev(:,1,:,:);

DGprev_west(:,2:end,:,:) = DGprev(:,1:(end-1),:,:) ;
DGprev_west(:,1,:,:) = DGprev(:,end,:,:);

DGprev_nort(:,:,1:(end-1),:) = DGprev(:,:,2:end,:) ;
DGprev_nort(:,:,end,:) = DGprev(:,:,1,:);

DGprev_sout(:,:,2:end,:) = DGprev(:,:,1:(end-1),:) ;
DGprev_sout(:,:,1,:) = DGprev(:,:,end,:);


DGprev_uppr(:,:,:,1:(end-1)) = DGprev(:,:,:,2:end);
DGprev_uppr(:,:,:,end) = DGprev(:,:,:,1);

DGprev_down(:,:,:,2:end) = DGprev(:,:,:,1:(end-1));
DGprev_down(:,:,:,1) = DGprev(:,:,:,end);



end


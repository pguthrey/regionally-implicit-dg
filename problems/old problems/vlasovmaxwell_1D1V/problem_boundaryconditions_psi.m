function [DGprev_east,DGprev_west,DGprev_nort,DGprev_sout] = problem_boundaryconditions_psi(DGprev,data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%load('data/data_parameters')

DGprev_east = NaN(data.thetaT,data.Nv1,data.Nv2);
DGprev_west = NaN(data.thetaT,data.Nv1,data.Nv2);
DGprev_nort = NaN(data.thetaT,data.Nv1,data.Nv2);
DGprev_sout = NaN(data.thetaT,data.Nv1,data.Nv2);

DGprev_east(:,1:(end-1),:) = DGprev(:,2:end,:);
    DGprev_east(:,end) = data.extrapolate_psi_x*DGprev(:,end);

DGprev_west(:,2:(end),:) = DGprev(:,1:(end-1),:);
    DGprev_west(:,1) = data.extrapolate_psi_x*DGprev(:,1);


end


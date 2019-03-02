function [ DGcoeffs_ghosts ] = problem_boundaryconditions_corrector(DGcoeffs,data)
%UNTITLED Summary of this function goes here

DGcoeffs_ghosts = zeros(data.thetaT,data.Nv1+2,data.Nv2+2,1);
DGcoeffs_ghosts(:,2:(data.Nv1+1),2:(data.Nv2+1),1) = DGcoeffs;

end


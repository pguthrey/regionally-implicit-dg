function [ DGcoeffs_ghosts ] = problem_boundaryconditions_psi(DGcoeffs,data)
%UNTITLED Summary of this function goes here

rx = data.rx_param;
ry = data.ry_param;

DGcoeffs_ghosts = zeros(data.thetaT,data.Nv1+2*rx,data.Nv2+2*ry,1);
DGcoeffs_ghosts(:,(rx+1):(data.Nv1+rx),(ry+1):(data.Nv2+ry),1) = DGcoeffs;

end


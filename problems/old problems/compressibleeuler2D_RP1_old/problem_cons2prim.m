function [ qprim ] = problem_cons2prim(qcons)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

qprim(4,:) = qcons(4,:);
qprim(3,:) = qcons(3,:)/qcons(1,:);
qprim(2,:) = qcons(2,:)/qcons(1,:);
qprim(1,:) = qcons(1,:);


end

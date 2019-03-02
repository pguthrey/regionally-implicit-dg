function [ qplot ] = problem_cons2plot(qcons,data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

qplot(6,:) = qcons(2,:);
qplot(5,:) = qcons(3,:)/qcons(1,:);
qplot(4,:) = (data.appdata.gamma-1)*(qcons(3,:)- qcons(2,:)^2/qcons(1,:)/2) ;
qplot(3,:) = qcons(3,:);
qplot(2,:) = qcons(2,:);
qplot(1,:) = qcons(1,:);


end


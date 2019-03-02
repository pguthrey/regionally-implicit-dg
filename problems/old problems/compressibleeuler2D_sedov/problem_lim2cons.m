function [qcons] = problem_lim2cons(qlim,data)

rho = qlim(1,:);
u = qlim(2,:);
v = qlim(3,:);
E = qlim(4,:);
p = qlim(5,:);
mx = qlim(6,:);
my = qlim(7,:); 
e = qlim(8,:);

qcons(4,:) = E;
qcons(3,:) = v*rho;
qcons(2,:) = u*rho;
qcons(1,:) = rho;


end


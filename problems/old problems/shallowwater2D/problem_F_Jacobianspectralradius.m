function [ rho ] = problem_F_Jacobianspectralradius(q,v1,v2)


%Shallow Water
alpha = 1/2;
h = q(1);
m = q(2);
p = q(3);
rho = max(abs(m/h+[-1,1]*sqrt(2*alpha*h)));
end


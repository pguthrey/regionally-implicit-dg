function [ rho ] = problem_F_Jacobianspectralradius(q,quadpoint,appdata)

u = q(2)/q(1);
N = u^2;
p = (appdata.gamma-1)*(q(3) - q(1)*N/2);
c = sqrt(appdata.gamma*p/q(1));

rho = abs(u)+c;

end


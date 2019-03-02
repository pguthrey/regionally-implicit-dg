function [ rho ] = problem_G_Jacobianspectralradius(q,quadpoint,appdata)

u = q(2)/q(1);
v = q(3)/q(1);
N = u^2+v^2;
p = (appdata.gamma-1)*(q(4) - q(1)*N/2);
c = sqrt(appdata.gamma*p/q(1));

rho = max(abs([v+c v-c]));

end


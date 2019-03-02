function [ drhodq ] = problem_F_drho_dq(q,quadpoint,appdata)

rho = q(1);
m = q(2);
E = q(3);
u = m/rho;
gamma = appdata.gamma;

dudq = [-m/rho^2 1/rho 0];
p = (gamma-1)*(E - m^2/rho/2);
dpdq = (gamma-1)*[u^2/2 -u 1];
dcdq = (gamma*p/rho).^(-1/2)/2.* ...
    gamma*(dpdq/rho - p/rho^2.*[1 0 0]);
    
drhodq = sign(u)*dudq + dcdq;

end


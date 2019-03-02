function [ drhodq ] = problem_F_drho_dq(q,quadpoint,appdata)

rho = q(1);
m = q(2);
E = q(3);
u = m/rho;
gamma = appdata.gamma;
N = u^2;
dudq = [-m/rho^2 1/rho 0];
p = (gamma-1)*(E - m^2/rho/2);

Z = p./rho;
dZdrho = (gamma-1)*(N/rho - E/rho^2);
dZdm = -(gamma-1)*m/rho^2;
dZdE = (gamma-1)/rho;

dcdq = sqrt(gamma)/sqrt(Z).*[dZdrho dZdm dZdE];
    
drhodq = sign(u)*dudq + dcdq;

end


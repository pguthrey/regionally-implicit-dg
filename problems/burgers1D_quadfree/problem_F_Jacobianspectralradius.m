function [ rho ] = problem_F_Jacobianspectralradius(q,~,appdata)

error('whoops?')
u = q(2)/q(1);
u2 = u*u;
p = (appdata.gamma-1)*(q(3) - q(1)*u2/2);
c = sqrt(appdata.gamma*p/q(1));

rho = abs(u)+c;

end


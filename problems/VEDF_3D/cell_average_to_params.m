function [params] = cell_average_to_params(n0,data)

Z = data.Z;
Ti = data.kbTi ;
Te = data.kbTe ;

%Constants
e = 1.60217662e-19; %Coulombs
epsilon_0 = 8.854187817620e-12;
m = 1.455e-25; %kg
gamma = 1;
xi = 0;

e0kTi = epsilon_0*data.kbTi;
e0kTe = epsilon_0*data.kbTe;

%derived parameters
ai = (0.75/pi/n0).^(1.0/3.0);
Gamma = e^2.0/(ai*e0kTi);
lambda_DE = sqrt(e0kTe/(4.0*pi*n0*(Z^2.0)*(e^2.0)));
kappa = ai/lambda_DE;
omega_p = sqrt(4.0*pi*e^2.0*n0/(m*epsilon_0));

%Lookup table for the excess energy
a = -0.899 + 0.5*kappa - 0.103*(kappa^2.0) + 0.003*(kappa^4.0);
b =  0.565             - 0.026*(kappa^2.0) - 0.003*(kappa^4.0);
c = -0.207             - 0.086*(kappa^2.0) + 0.018*(kappa^4.0);
d = -0.031             + 0.042*(kappa^2.0) - 0.008*(kappa^4.0);
Efit = a*Gamma+b*(Gamma^(1.0/3.0))+c+d*(Gamma^(-1.0/3.0));

%
eta_kappa_raw = [	0.1	,	0.1	,	0.1	,	0.1	,	0.1	,	0.1	,	0.5	,	0.5	,	0.5	,	0.5	,	0.5	,	0.5	,	1	,	1	,	1	,	1	,	1	,	1	];
eta_Gamma_raw =	[	2	,	5.02	,	10	,	20	,	50	,	100	,	2	,	5.01	,	10	,	19.9	,	50.2	,	100	,	2	,	4.99	,	9.9	,	19.8	,	49.4	,	99	];
etastar_raw =	[	0.502	,	0.128	,	0.0686	,	0.0691	,	0.0912	,	0.206	,	0.5	,	0.13	,	0.0874	,	0.0629	,	0.0861	,	0.192	,	0.486	,	0.172	,	0.106	,	0.0904	,	0.0964	,	0.179	];

etastar = griddata(eta_kappa_raw,eta_Gamma_raw,etastar_raw,kappa,Gamma);



eta = etastar*m*n0*omega_p*ai^2.0;
etatilde = (4/3*eta+xi)/(m*n0*omega_p*ai^2);

dEfit_dgammatilde = a+b/3.0*Gamma^(-2.0/3.0)-d/3.0*Gamma^(-4.0/3.0);
Ectilde = Efit; %unitless Ec?  *kTi
Ec = Ectilde*e0kTi;
dEc_dGammatilde = dEfit_dgammatilde - Efit/(3.0*Gamma);
mu = 1.0+Ectilde/3.0+Gamma/9.0*dEc_dGammatilde;
psi = - Ectilde/15.0 - Gamma/9.0*dEc_dGammatilde;

tau = (3*etatilde*Gamma)/psi/omega_p;

    
params.Gamma = Gamma;
params.kappa = kappa;
params.tau = tau;
params.xi = xi;
params.mu = mu;
params.eta = eta;


end


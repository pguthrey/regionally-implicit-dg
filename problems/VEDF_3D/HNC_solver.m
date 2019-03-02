clc
clear all
close all

Gamma = 2
kappa = .1

%N_bins = 128
N_bins = 512

%--------------------------------------------------------------------------------------------
% run parameters

num_iterations = 1000
epsilon = 0.95
alpha = 1.0

% r binning
R_max = 25.0
del_r = R_max/N_bins
% error is minimized when using bin centers
r_array = linspace(del_r/2,R_max-del_r/2,N_bins)
energy_iter = zeros(num_iterations,1)
pressure_iter = zeros(num_iterations,1)

% k binning
del_k = pi/((N_bins+1)*del_r)
K_max = del_k*N_bins
k_array = linspace(del_k/2,K_max-del_k/2,N_bins)


% precompute factors
fact_r_2_k = 2*pi*del_r
fact_k_2_r = del_k/(4.*pi^2)
dimless_dens = 3./(4*pi)

%--------------------------------------------------------------------------------------------

wash = @(X) double(py.array.array('d',py.numpy.nditer(X)))

% f(k) = FT[f(r)]
FT_r_2_k = @(input_array) fact_r_2_k.*wash(py.scipy.fftpack.dst(r_array.*input_array))./k_array;

% f(r) = IFT[f(k)]
FT_k_2_r = @(input_array) fact_k_2_r.*wash(py.scipy.fftpack.idst(k_array.*input_array))./r_array;

%--------------------------------------------------------------------------------------------
%Gamma = params.Gamma
%kappa = params.kappa
%--------------------------------------------------------------------------------------------

% Returns a tuple of arrays with [u(r),u_s(r)]; it is the short-ranged u_s(r) that is used in the HNC solver.
% Note that u(r) itself is not needed for HNC, but is used for computing the energy, which is used to test convergence.
% This function assumes global variables Gamma, kappa and alpha.
u1 = @(r_in) Gamma*exp(-kappa*r_in)./r_in ; 
u2 = @(r_in) Gamma*exp(-(alpha + kappa)*r_in)./r_in ;

% when du/dr is known, precompute it here; only needed for computing the pressure, not for the HNC solution itself
deriv_u_r = @(r_in)  -Gamma*exp(-kappa*r_in).*(1 + kappa*r_in)./r_in.^2;

% The initial condition is given through c(k), defined here. This should be consistent with the other definitions
% in this cell.
initial_c_k = @(k_in) -4*pi*Gamma./(k_array.^2 + kappa^2);

% Long range part of the potential in Fourier space. This could be computed with the Fourier routines, but
% I chose to precopmute things analytically when possible.
u_long_l = @(k_in) 4*pi*Gamma*(alpha^2 + 2*alpha*kappa)./((k_array.^2 + kappa^2).*(k_array.^2+(alpha + kappa)^2));

%--------------------------------------------------------------------------------------------
u_r = u1(r_array);
u_s_r = u2(r_array);
d_u_d_r = deriv_u_r(r_array);
c_k = initial_c_k(k_array);
u_l_k = u_long_l(k_array);
c_s_k = c_k + u_l_k;
S_k = 1+dimless_dens*c_k./(1-dimless_dens*c_k);
% start the iteration loop
for iteration = 1:num_iterations
    gamma_s_k = (dimless_dens*c_s_k.*c_k - u_l_k)./(1-dimless_dens*c_k);
    gamma_s_r = FT_k_2_r(gamma_s_k);
    g_r = exp(gamma_s_r - u_s_r);
    
    energy_iter(iteration) = 1.5*del_r*sum(r_array.^2.*u_r.*(g_r-1));
    pressure_iter(iteration) = -0.5*del_r*sum(r_array.^3.*d_u_d_r.*(g_r-1));
    new_c_s_r = g_r - 1 - gamma_s_r;
    c_s_k = FT_r_2_k(new_c_s_r);
    c_k = c_s_k - u_l_k;
end

g_data = [r_array;g_r]';  % g(r)
c_data = [k_array;c_k]';  % c(k)
S_data = [k_array;S_k]';  % S(k)

%Let C(|x-y|) = c(|x-y|)|x-y|  
% C(x)        = c(|x|)|x|
% dFdn        = int n(x)C(|x-y|) dy 
% F[dFdn](k)  = F[n](k)F[C](k)    
% dFdn(x)     = F^-1[ F[n](k)F[C](k) ](x)
c_x = k_array;         % x=k values
c_f = FT_k_2_r(c_k);   % c real space
c_F = FT_r_2_k(c_f);  % c in fourier space
c_fr = c_f.*c_x;        % c*r in real space
c_Fr = FT_r_2_k(c_fr); % c*r in fourier space


figure
plot(c_x,c_f)
figure
plot(g_data(:,1),g_data(:,2))
figure
plot(S_data(:,1),S_data(:,2))
%return c_x,c_f,c_F,c_fr,c_Fr,energy_iter(end),pressure_iter(end)

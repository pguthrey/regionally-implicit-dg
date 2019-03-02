clc
clear all
close all

addpath(['problems/' problemname])
addpath('testfunctions')
addpath('initialize')


M = 4%[4 6]
r = [1]
data.space_dims = 1;

predbasispick = 1%[1 2]
corrbasispick = 2%[2 1]
data.basis = 'canonical';
data.M = M;
%data.predictor_solver = 'NewtonIteration';
switch predbasispick
    case 1
        data.predictorbasis = 'Q';
    case 2
        data.predictorbasis = 'P';
end

switch corrbasispick
    case 1
        data.correctorbasis = 'Q';
    case 2
        data.correctorbasis = 'P';
end

data.methodtype = 'implicit';
%[data] = initialize_method(data);
[data] = initialize_Basis(data,data.space_dims);
%[data] = initialize_gauss_quadrature(data,data.M);
%[data] = initialize_innerproduct_quadrature(data);
%[data] = initialize_mesh(data);

syms t real
syms x real
vectpsi = sym(zeros(thetaT,Neqns));
vectphi = sym(zeros(theta,Neqns));
for kay = 1:thetaT
    for eqn = 1:Neqns
        ell = sysBO(eqn,kay);
        coeffs = sysBOcoeffs(eqn,kay);
        vectpsi(kay,eqn) = coeffs.*testfunction_psi([t,x,NaN,NaN],ell,data);
    end
end
for kay = 1:theta
    for eqn = 1:Neqns
        ell = sysBOs(eqn,kay);
        coeffs = sysBOscoeffs(eqn,kay);
        vectphi(kay,eqn) = testfunction_phi([x,NaN,NaN],ell,data);
    end
end



Ue = sym('u_east',[data.thetaTtilde,1]);
Us = sym('u_star',[data.thetaTtilde,1]);
Uw = sym('u_west',[data.thetaTtilde,1]);
u_east = vectpsi'*Ue;
u_star = vectpsi'*Us;
u_west = vectpsi'*Uw;
 
for i = thetaT:-1:1
        str = ['u_star' num2str(i) '= qstar(' num2str(i) ');'];
        disp(str)   
end
%{
for i = thetaT:-1:1
        str = ['u_west' num2str(i) '= qwest(' num2str(i) ');'];
        disp(str)   
end
for i = thetaT:-1:1
        str = ['u_east' num2str(i) '= qeast(' num2str(i) ');'];
        disp(str)   
end
%}
 
Flux = int(subs(vectpsi,x,1)*subs(u_star^2/2,x,1),t,-1,1);
thing = 'trunceast_temp';
for i = thetaT:-1:1
    if ~(Flux(i) == 0)           
        str = [thing '_pre(' num2str(i) ',1) = ' char(vpa(Flux(i),16)) ';'];
        disp(str)            
    end
end
str = [thing ' = + nuv1*' thing '_pre;'];
disp(str)  
 
Flux = int(subs(vectpsi,x,-1)*subs(u_star^2/2,x,-1),t,-1,1);
thing = 'truncwest_temp';
for i = thetaT:-1:1
    if ~(Flux(i) == 0)           
        str = [thing '_pre(' num2str(i) ',1) = ' char(vpa(Flux(i),16)) ';'];
        disp(str)            
    end
end
str = [thing ' = - nuv1*' thing '_pre;'];
disp(str)
 
F_dql = subs(u_star,x,1) + (subs(u_star,x,1)+subs(u_east,x,-1))/2 ;
Jacobian = int(subs(vectpsi,x,1)*F_dql*subs(vectpsi',x,1),t,-1,1);
thing = 'Jac_east_cell_temp';
for i = thetaT:-1:1
    for j = thetaT:-1:1
        if ~(Jacobian(i,j) == 0)           
            str = [thing '_pre(' num2str(i) ',' num2str(j) ') = ' char(vpa(Jacobian(i,j),16)) ';'];
            disp(str)       
        end
    end
end
str = [thing ' = nuv1*' thing '_pre;'];
disp(str)
 
Jacobian = int(subs(vectpsi*u_star,x,1)*subs(vectpsi',x,1),t,-1,1);
thing = 'Jac_trunc_east_temp';
for i = thetaT:-1:1
    for j = thetaT:-1:1
        if ~(Jacobian(i,j) == 0)           
            str = [thing '_pre(' num2str(i) ',' num2str(j) ') = ' char(vpa(Jacobian(i,j),16)) ';'];
            disp(str)       
        end
    end
end
str = [thing ' = nuv1*' thing '_pre;'];
disp(str)
 
Jacobian = int(subs(vectpsi*u_star,x,-1)*subs(vectpsi',x,-1),t,-1,1);
thing = 'Jac_trunc_west_temp';
for i = thetaT:-1:1
    for j = thetaT:-1:1
        if ~(Jacobian(i,j) == 0)           
            str = [thing '_pre(' num2str(i) ',' num2str(j) ') = ' char(vpa(Jacobian(i,j),16)) ';'];
            disp(str)
        end
    end
end
str = [thing ' = -nuv1*' thing '_pre;'];
disp(str)
 
Flux = int(int(diff(vectpsi,x)*u_star^2/2,x,-1,1),t,-1,1);
thing = 'residual_cell';
for i = thetaT:-1:1
    if ~(Flux(i) == 0)           
        str = [thing '_pre(' num2str(i) ',1) = ' char(vpa(Flux(i),16)) ';'];
        disp(str)            
    end
end
str = 'residual_cell = residual_cell - nuv1*residual_cell_pre;';
disp(str)
 
Jacobian = int(int(diff(vectpsi,x)*u_star*vectpsi',x,-1,1),t,-1,1);
for i = thetaT:-1:1
    for j = thetaT:-1:1
        if ~(Jacobian(i,j) == 0)           
            str = ['Jac_cell_pre(' num2str(i) ',' num2str(j) ') = ' char(vpa(Jacobian(i,j),16)) ';'];
            disp(str)            
        end
    end
end
str = 'Jac_cell = data.Psi_futr - nuv1*Jac_cell_pre;';
disp(str)
 
     
keyboard
%}
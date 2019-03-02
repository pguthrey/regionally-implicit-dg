clc
clear all
close all

dims = 1

% 1D 
syms n m1  Mx
fx = [ m1	;  m1^2/n; 	m1/n*Mx ];
vars = [n ; m1 ; Mx ];
for j = 1:3   
    Ax(:,j) = diff(fx,vars(j));
end

syms dpdn
Ax(2,1) = Ax(2,1) + dpdn;

disp('Ax eigs')
[V,D] = eig(Ax)


return

disp('--------------------------------------------------------------')

% 2D
clear Ax 
syms  m2   My  
fx = [ m1	;  m1^2/n; 	m1*m2/n;	m1/n*Mx ;	m1/n*My ; ];
fy = [ m2	;  m1*m2/n; m2*m2/n;    m2/n*Mx ;	m2/n*My ; ];
vars = [n ; m1 ; m2 ; Mx ; My ];
for j = 1:5    
    Ax(:,j) = diff(fx,vars(j));
    Ay(:,j) = diff(fy,vars(j));
end

syms dpdn
Ax(2,1) = Ax(2,1) + dpdn;
Ay(3,1) = Ay(3,1) + dpdn;

disp('Ax eigs')
[V,D] = eig(Ax)

disp('Ay eigs')
[V,D] = eig(Ay)


disp('--------------------------------------------------------------')


% 3D
clear Ax  Ay
syms  m3  Mz 
fx = [ m1	;  m1^2/n; 	m1*m2/n;	m1*m3/n ;	 m1/n*Mx ;	m1/n*My ;	m1/n*Mz ];
fy = [ m2	;  m1*m2/n; 	m2*m2/n;	m2*m3/n ;	 m2/n*Mx ;	m2/n*My ;	m2/n*Mz ];
fz = [ m3	;  m1*m3/n; 	m2*m3/n;	m3*m3/n ;	 m3/n*Mx ;	m3/n*My ;	m3/n*Mz ];
vars = [n ; m1 ; m2 ; m3 ; Mx ; My ; Mz ];
for j = 1:7    
    Ax(:,j) = diff(fx,vars(j));
    Ay(:,j) = diff(fy,vars(j));
    Az(:,j) = diff(fz,vars(j));
end

syms dpdn
Ax(2,1) = Ax(2,1) + dpdn;
Ay(3,1) = Ay(3,1) + dpdn;
Az(4,1) = Az(4,1) + dpdn;

disp('Ax eigs')
[V,D] = eig(Ax)

disp('Ay eigs')
[V,D] = eig(Ay)


disp('Az eigs')
[V,D] = eig(Az)




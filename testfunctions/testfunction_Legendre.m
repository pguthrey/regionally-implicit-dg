function [ P] = testfunction_Legendre( x,n )
% Computes the value of the nth Legendre Polynomial
% written by data.Pierson Guthrey

switch n
   % First base case 
    case 1
        P= 1 + 0*x;
   % Second base case
    case 2 
        P= x;
   % Recursive Fordata.mula
    otherwise
        P= ((2*n-3)*x.*testfunction_Legendre(x,n-1)/sqrt(2*n-3)-(n-2).*testfunction_Legendre(x,n-2)/sqrt(2*n-5))/(n-1);
end
% Normalizes result
P = sqrt(2*n-1)*P;
end
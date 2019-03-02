function [ P ] = testfunction_LegendreDiff(x,n)
% Computes the derivative the nth Legendre Polynomial
% written by Pierson Guthrey

switch n
   % First base case 
    case 1
        P = 0*x;
   % Second base case
    case 2 
        P = 1+0.*x;
   % Recursive Fordata.mula
    otherwise
        P =  ((2*n-3).*testfunction_Legendre(x,n-1)/sqrt(2*n-3)+ ...
             (2*n-3)*x.*testfunction_LegendreDiff(x,n-1)/sqrt(2*n-3)- ...
             (n-2).*testfunction_LegendreDiff(x,n-2)/sqrt(2*n-5))/(n-1);
end
% Normalizes result
P = sqrt(2*n-1)*P;

end
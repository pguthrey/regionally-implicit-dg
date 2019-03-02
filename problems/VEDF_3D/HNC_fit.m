function [coeff,scaling] = HNC_fit(xvals,cvals)

capprox = @(r_array,params) params(1)*exp(-params(2)*xvals.^2.0);

del_r = xvals(2)-xvals(1);
objective = @(params) norm(capprox(xvals,params) - cvals,2)*sqrt(del_r);

A_guess = -10 ;
a_guess = .8  ;
            
guess = [A_guess,a_guess]       ;

options = optimoptions('fminunc','Display','off');

soln = fminunc(objective,guess,options) ;
          
coeff = soln(1);
scaling = soln(2);

end


function [searchdir] = predictor_backslash_solver(Jacobian,residual,thetaT)
     
    coeffs = (1:thetaT)+(5-1)*thetaT;
    soln = Jacobian\residual;
    searchdir = soln(coeffs);

    
end


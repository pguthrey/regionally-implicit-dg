function [ ] = problem_plotting_IC(solution,data)
% Computes the exact solution
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    (t,v1,v2,v3) : location in space-time 
% OUTdata.PUTS   exact   : the given exact solution evaluated at the location
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

plot_1D_system(solution,0,data)

end
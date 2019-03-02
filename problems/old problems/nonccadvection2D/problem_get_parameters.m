function [data] = problem_get_parameters(data)

data.space_dims = 2;

%%{
data.Fspeedmax = 1;
data.Gspeedmax = 1;
data.Hspeedmax = 0;
%}
data.measureerror = true;
data.boundsv1v2v3 = [-5 5 -5 5 0 0];

data.linearcase = 'linear';
 
data.Tfinal = 2*pi;
%NonCCAdvection
data.boundsq  = [-.05 1.05];
data.Neqns = 1;
data.varname{1} = 'Mass';
data.Nplotvars = 1;
data.plotname{1} = 'Mass';
data.boundsplot  = [-.05 1.05];

data.space_dims = 2;
data.appdata = 0;
end


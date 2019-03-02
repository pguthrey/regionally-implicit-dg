function [data] = problem_get_parameters(data)

data.space_dims = 2;

data.boundsv1v2v3 = [0 2*pi 0 2*pi 0 0];

data.Tfinal = .4;

%NonCCAdvection
data.boundsq  = [-.1 1.1];
data.Neqns = 1;
data.varname{1} = 'Velocity';
%data.formulation = 'strong';

data.Fspeedmax = 1.05;
data.Gspeedmax = 1.05;
data.Hspeedmax = 0;

data.appdata = 0;

data.Nplotvars = 1;
data.plotname{1} = 'Velocity';

data.boundsplot  = [-.1 1.1];
             
data.linearcase = 'nonlinear';

data.measureerror = true;

end


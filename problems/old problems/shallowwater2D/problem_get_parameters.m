function [data] = problem_get_parameters(data)

data.space_dims = 2;

data.measureerror = true;
data.boundsv1v2v3 = [-1 1 -1 1 0 0];

error('shock time?')
data.Tfinal = pi/7;

%NonCCAdvection
data.boundsq  = [-.1 1.1 -1.1 1.1 -1.1 1.1];
data.Neqns = 1;
data.varname{1} = 'hieght';
data.varname{2} = 'x-velocity';
data.varname{3} = 'y-velocity';
%data.formulation = 'strong';

data.Fspeedmax = 1;
data.Gspeedmax = 1;
data.Hspeedmax = 1;

end


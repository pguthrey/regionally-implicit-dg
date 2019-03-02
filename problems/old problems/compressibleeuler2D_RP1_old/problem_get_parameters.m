function [data] = problem_get_parameters(data)

data.space_dims = 2;

data.measureerror = false;
data.boundsv1v2v3 = [0 1 0 1 0 0];

data.Tfinal = .8;

%NonCCAdvection
data.boundsq  = [-.1 2
                 -2 2
                 -2 2
                 -.1 3];
data.Neqns = 4;
data.varname{1} = 'Mass';
data.varname{2} = 'x momentum ';
data.varname{3} = 'y momentum';
data.varname{4} = 'Energy';
%data.formulation = 'strong';

data.appdata.gamma = 5/3; 

p = 1;
c = 1;

data.Fspeedmax = 3;
data.Gspeedmax = 3;
data.Hspeedmax = 0;

end


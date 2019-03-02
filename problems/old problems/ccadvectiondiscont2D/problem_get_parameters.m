function [data] = problem_get_parameters(data)

data.space_dims = 2;

data.measureerror = true;
data.boundsv1v2v3 = [-1 1 -1 1 0 0];

data.Tfinal = pi/7;

%NonCCAdvection
data.boundsq  = [-.02 1.02];
data.Neqns = 1;
data.varname{1} = 'Mass';
%data.formulation = 'strong';

data.Nplotvars = 1;
data.plotname{1} = 'Mass';

data.boundsplot  = [-.1 1.1 ];
             
             
data.appdata.nuf = 0;
data.appdata.nug = 1;
data.appdata.nuh = 0;

data.Fspeedmax = abs(data.appdata.nuf);
data.Gspeedmax = abs(data.appdata.nug);
data.Hspeedmax = abs(data.appdata.nuh);

end


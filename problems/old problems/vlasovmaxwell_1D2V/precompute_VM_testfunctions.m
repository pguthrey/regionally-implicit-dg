function [data] = precompute_VM_testfunctions(data)
% data.Precomputes the test functions at quadrature points
% written by data.Pierson Guthrey
% -------------------------------------------------
% INdata.PUTS    
% OUTdata.PUTS   data.Precompute test funtions
% Note: other variables may be %loaded in from the problem files
% ------------------------------------------------------------------------

%load('vlasovmaxwell/data.appdata')

zeta1D2V = zeros(1,data.theta1D2V,data.P,data.P,data.P);
zeta1D = zeros(1,data.theta1D,data.P);
zeta1Dplot = zeros(1,data.theta1D,data.P);

data.vmplotlocs = linspace(-1,1,data.P+2);

for kay = 1:data.theta1D
    for xquad = 1:data.P
        xii = data.locs(xquad);
        data.zeta1D(1,kay,xquad) = testfunction_zeta1D(xii,kay,data);
    end
end

for kay = 1:data.theta1D
    for xquad = 1:data.P+2
        xii = data.vmplotlocs(xquad);
        data.zeta1Dplot(1,kay,xquad) = testfunction_zeta1D(xii,kay,data);
    end
end

for kay = 1:data.theta1D2V  
    for xquad = 1:data.P+2
    for v1quad = 1:data.P+2
        xii = data.vmplotlocs(xquad);
        eta1 = data.vmplotlocs(v1quad);
        data.zeta1D1Vplot(1,kay,xquad,v1quad) = testfunction_zeta1D2V(xii,eta1,0,kay,data);
    end
    end
end

for kay = 1:data.theta1D2V  
    for xquad = 1:data.P
    for v1quad = 1:data.P
    for v2quad = 1:data.P
        xii = data.locs(xquad);
        eta1 = data.locs(v1quad);
        eta2 = data.locs(v2quad);
        data.zeta1D2V(1,kay,xquad,v1quad,v2quad) = testfunction_zeta1D2V(xii,eta1,eta2,kay,data);
    end
    end
    end
end
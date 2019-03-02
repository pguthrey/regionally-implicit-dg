function [data] = initialize_innerproduct_quadrature(data)
% Initialized the inner product quadrature 
% written by Pierson Guthrey

D = data.space_dims;
P = data.P;
data.Pplot = data.M+1;
data.Plimiterpsi = data.P+2;
data.Plimiter = data.P+2;

data.plotlocs = linspace(-1,1,data.Pplot);
data.limiterlocs = [-1;data.locs;1];

data.Pdm1 = P^(D-1);
data.Pd = P^(D);
data.Pdp1 = P^(D+1); 
data.Pdplot = data.Pplot^D; 
data.Pdlimiter = data.Plimiter^D; 
data.Pdplotpsi = data.Pplot^(D+1); 
data.Pdlimiterpsi = data.Plimiterpsi^(D+1); 

Dm1base =(1:P^(D-1))';
Dbase =(1:P^(D))';
Dp1base =(1:P^(D+1))';
plotbase = (1:data.Pplot^D)';
limiterbase = (1:data.Plimiter^D)';
limiterpsibase = (1:data.Plimiterpsi^D)';
plotbasepsi = (1:data.Pplot^(D+1))';

data.Dm1list = zeros(P^(D-1),max([D-1 1]));
data.Dlist = zeros(P^D,D); 
data.Dp1list = zeros(P^(D+1),D+1);
data.plotlist = zeros(data.Pplot^D,D);
data.limiterlist = zeros(data.Plimiter^D,D);
data.psiplotlist = zeros(data.Pplot^(D+1),D+1);
data.psilimiterlist = zeros(data.Plimiterpsi^(D+1),D+1);

%data.plotlistpsi = zeros(P^(D+1),D+1);

Dp1quadwgts = zeros(data.Pdp1,D+1);
Dm1quadwgts = zeros(data.Pdm1,max([D-1 1]));
Dquadwgts = zeros(data.Pd,D);

for j=1:max([D-1 1])
    data.Dm1list(:,j)=mod(ceil(Dm1base./P.^(j-1))-1,P)+1;
end
for j=1:D
    data.Dlist(:,j)=mod(ceil(Dbase./P.^(j-1))-1,P)+1;
end
for j=1:(D+1)
    data.Dp1list(:,j)=mod(ceil(Dp1base./P.^(j-1))-1,P)+1;
    data.plotlistpsi(:,j)=mod(ceil(plotbasepsi./(data.Pplot).^(j-1))-1,data.Pplot)+1;
    data.limiterlistpsi(:,j)=mod(ceil(limiterpsibase./(data.Plimiterpsi).^(j-1))-1,data.Plimiterpsi)+1;
end

for j=1:D
    data.plotlist(:,j)= mod(ceil(plotbase./(data.Pplot).^(j-1))-1,data.Pplot)+1;
    data.limiterlist(:,j)= mod(ceil(limiterbase./(data.Plimiter).^(j-1))-1,data.Plimiter)+1;
end

if data.space_dims >= 2    
    for d = 1:max([D-1 1])
    for k = 1:data.Pdm1
        data.Dm1quadlocs(k,d) = data.locs(data.Dm1list(k,d));
        Dm1quadwgts(k,d) = data.wgts1D(data.Dm1list(k,d));
    end
    end
else
    data.Dm1quadlocs = 0;
    Dm1quadwgts = 1;
end

for d = 1:D
for k = 1:data.Pd
    data.Dquadlocs(k,d) = data.locs(data.Dlist(k,d));
    Dquadwgts(k,d) = data.wgts1D(data.Dlist(k,d));
end
end

for d = 1:(D+1)
for k = 1:data.Pdp1
    data.Dp1quadlocs(k,d) = data.locs(data.Dp1list(k,d));
    Dp1quadwgts(k,d) = data.wgts1D(data.Dp1list(k,d));
end
end

for d = 1:D
for k = 1:data.Pdplot
    data.Dplotlocs(k,d) = data.plotlocs(data.plotlist(k,d));
end
end

for d = 1:D
for k = 1:data.Pdlimiter
    data.Dlimiterlocs(k,d) = data.limiterlocs(data.limiterlist(k,d));
end
end

for d = 1:(D+1)
for k = 1:data.Pdplotpsi
    data.Dplotlocspsi(k,d) = data.plotlocs(data.plotlistpsi(k,d));
end
end


data.Dm1quadwgts = prod(Dm1quadwgts,2);
data.Dquadwgts = prod(Dquadwgts,2);
data.Dp1quadwgts = prod(Dp1quadwgts,2);

end
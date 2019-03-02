function [ temp ] = integrate_fnautds(sigma,t,x,v1,v2,data)
%{
temp = 0;
for tquad = 1:data.P    
    tau = (data.locs(tquad)+1)/2*t;
    wgt = data.wgts1D(tquad);    
    temp = temp + sigma(tau,x,v1,v2)*wgt;
end
%}

%Trapezoid rule
temp = (sigma(0,x,v1,v2) + sigma(t,x,v1,v2))/2;




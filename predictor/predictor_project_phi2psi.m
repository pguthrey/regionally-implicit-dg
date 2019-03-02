function [ DGpreprediction ] = predictor_project_phi2psi(DGprevioussolution,data)
% written by Pierson Guthrey

thetaT = data.thetaT;
Nv1 = data.Nv1;
Nv2 = data.Nv2;
Nv3 = data.Nv3;
DGpreprediction  = zeros(thetaT,Nv1,Nv2,Nv3);
projection = data.PHI2PSI;

%parfor iv1 = 1:Nv1
for iv1 = 1:Nv1
    DGpreprediction_temp = NaN(thetaT,Nv2,Nv3);    
for iv2 = 1:Nv2
for iv3 = 1:Nv3
    phicoeffs = DGprevioussolution(:,iv1,iv2,iv3);
    DGpreprediction_temp(:,iv2,iv3) = projection*phicoeffs;
end
end     
    DGpreprediction(:,iv1,:,:) = DGpreprediction_temp;
end

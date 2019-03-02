function [data] = initialize_Basis(data,D)
% Creates the matrices that contain information about the DG basis
% written by Pierson Guthrey
% -------------------------------------------------
% INPUTS    data.M
%           D
%           basis       : choice of basis type
% OUTPUTS   data.BO          : matriv1 of 2D spacetime basis orders
%                           rows are basis function orders, 
%                           1st col is order of t basis, 
%                           2,3,... are orders of spatial basis 
%           data.BOs         : matriv1 of 2D space basis orders
%                           rows are basis function orders,
%                           1,2,... cols are orders of spatial basis       
%           L           : data.number of elements in the spacetime basis
%           Ls          : data.number of elements in the space basis
%           data.sysBO       : matriv1 of 2D spacetime vector basis orders
%           data.sysBOs      : matriv1 of 2D space vector basis orders
%           data.thetaT       : data.number of elements in the spacetime vector basis
%           data.theta      : data.number of elements in the space vector basis
%           data.sysBOcoeffs : coefficients of the 2D spacetime vector basis
%           data.sysBOscoeffs: coefficients of the 2D space vector basis
% Note: other variables may be %loaded in from the problem parameter files
% ------------------------------------------------------------------------

J=(1:data.M^(D+1))';
LM=zeros(data.M^(D+1),D+1); %cartesian product of bases
for j=1:D+1
    LM(:,D+1-j+1)=mod(ceil(J./data.M.^(j-1)),data.M);
end
LM(LM==0)=LM(LM==0)+data.M;
LO= LM-1; %Actual polynomial orders

%%{ 
%old way of doing this
switch data.predictorbasis
    case 'P'
        data.BO = LM(sum(LO,2)<=data.M-1,:);%Basis orders picked from polynomial orders <= data.M-1
        L=length(data.BO(:,1)); 
    case 'Q'
        data.BO = LM;%(max(LO,[],2)<=data.M-1,:);%Basis orders picked from polynomial orders <= data.M-1
        L=length(data.BO(:,1)); 
end

BOs_all = LM(LM(:,1)==1,2:end);
switch data.correctorbasis
    case 'P'                
        BOs = BOs_all(sum(BOs_all-1,2)<=data.M-1,:);
        Ls = length(BOs(:,1));
        data.BOs= BOs; 
    case 'Q'
        BOs = BOs_all;
        Ls = length(BOs(:,1));
        data.BOs= BOs;%(1:Ls,2:end);
end
%}

clear BOs

switch data.space_dims
    case 1
        cnt = 0;
        num_basis_cmpts_1d = data.M;
        for order = 1:num_basis_cmpts_1d
            BOs(order,:) = [order];
        end
    case 3
        cnt = 0;
        num_basis_cmpts_1d = data.M;
        for order = 1:num_basis_cmpts_1d
            for i = 1:num_basis_cmpts_1d
                for j = 1:num_basis_cmpts_1d
                    for k = 1:num_basis_cmpts_1d
                        if ((i+j+k-2)==order)
                            cnt = cnt+1;
                            BOs(cnt,:) = [i j k];
                        end
                    end
                end
            end
        end
end

%data.BOs= BOs;

% New way of doing this
%{
switch data.predictorbasis
    case 'P'        
        data.BO = LM(sum(LM == data.M,2) <= 1,:);%Basis orders picked from polynomial orders <= data.M-1
        %data.BO = LM(sum(LO,2)<=data.M-1,:);%Basis orders picked from polynomial orders <= data.M-1
        L=length(data.BO(:,1)); 
    case 'Q'
        data.BO = LM;%(max(LO,[],2)<=data.M-1,:);%Basis orders picked from polynomial orders <= data.M-1
        L=length(data.BO(:,1)); 
end

BOs_all = LM(LM(:,1)==1,2:end);
switch data.correctorbasis
    case 'P'                
        %BOs = BOs_all(sum(BOs_all-1,2)<=data.M-1,:);
        BOs = BOs_all(sum(BOs_all == data.M,2) <= 1,:);
        Ls = sum(BOs(:,1)==1);
        data.BOs= BOs; 
    case 'Q'
        BOs = BOs_all;
        Ls = length(BOs(:,1));
        data.BOs= BOs;%(1:Ls,2:end);
end
keyboard
%}


data.L = L;
data.Ls = Ls;
data.thetaT = L*data.Neqns;
data.theta = Ls*data.Neqns;
data.sysBO = zeros(data.Neqns,L*data.Neqns);
data.sysBOcoeffs = zeros(data.Neqns,L*data.Neqns);
data.sysBOs = zeros(data.Neqns,Ls*data.Neqns);
data.sysBOscoeffs = zeros(data.Neqns,Ls*data.Neqns);

if strcmp(data.basis,'canonical')
    for i = 1: data.Neqns
        J = (i-1)*L+(1:L);
        data.sysBO(i,J) = (1:L);         
        data.sysBOcoeffs(i,J) = ones(1,L);         
        J = (i-1)*Ls+(1:Ls);
        data.sysBOs(i,J) = (1:Ls);         
        data.sysBOscoeffs(i,J) = ones(1,Ls);         
    end
else
    disp('Invalid basis choice')
end


end



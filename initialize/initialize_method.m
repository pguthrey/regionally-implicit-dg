function [data] = initialize_method(data)
% written by Pierson Guthrey

if data.r_param == 0
    data.method = 'LIDG';
else
data.method = ['RI' num2str(data.r_param) 'DG'];
end

data.rx_param = data.r_param;
data.ry_param = 0;
data.rz_param = 0;

if data.space_dims >= 2
    data.ry_param = data.r_param;
    if data.space_dims >= 3
        data.rz_param = data.r_param;
    end
end

data.region_per_dimension = 2*data.r_param +1;
data.cells_per_region = (data.region_per_dimension)^data.space_dims;
data.main_cell = (data.cells_per_region+1)/2;

switch data.methodtype
    case 'implicit'
        [data] = initialize_method_implicit(data);
end

if data.cfl == 0
   warning('CFL not selected!, using pessimistic CFL') 
   data.cfl = 0.05;
end


end       
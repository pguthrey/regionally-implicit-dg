function [data] = initialize_gauss_hermite_quadrature(data,order)
% written by Pierson Guthrey

switch order
    case 1
        mat = [1	0.0         1.77245	1.77245];
    case 2
        mat = [1	-0.707107	0.886227	1.46114
                2	0.707107	0.886227	1.46114];
    case 3
        mat = [ 1	-1.22474	0.295409	1.32393
                2	0.0         1.18164     1.18164
                3	1.22474     0.295409	1.32393];
    case 4
        mat= [  1	-1.65068	0.0813128	1.24023
                2	-0.524648	0.804914	1.05996
                3	0.524648	0.804914	1.05996
                4	1.65068     0.0813128	1.24023];         
       otherwise
            error('Gauss Hermite quadrature not implemented')
end

data.GH_wgts1D = mat(:,2); 
data.GH_locs = mat(:,4); 


end


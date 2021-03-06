function [data] = initialize_gauss_quadrature(data,order)
% written by Pierson Guthrey


switch order
    case 1
         wgts = 2;
         data.locs = 0;
    case 2
         wgts = [    1
                    1]; 
         data.locs = [  -0.5773502691896257
                    0.5773502691896257];
    case 3
         wgts = [  0.5555555555555556 
                  0.8888888888888888
                  0.5555555555555556 ];
         data.locs = [  -0.7745966692414834 
                   0
                   0.7745966692414834];
    case 4
         wgts = [ 0.3478548451374538
                 0.6521451548625461                 
                 0.6521451548625461
                 0.3478548451374538];
         data.locs = [ -0.8611363115940526
                -0.3399810435848563                
                0.3399810435848563
                0.8611363115940526];
    case 5
         wgts = [0.2369268850561891
                0.4786286704993665	
                0.5688888888888889
                0.4786286704993665
                0.2369268850561891];
         data.locs = [-0.9061798459386640
                -0.5384693101056831
                0
                0.5384693101056831
                0.9061798459386640];
    case 6
           wgts = [0.1713244923791704
                0.3607615730481386	
                0.4679139345726910
                0.4679139345726910
                0.3607615730481386	
                0.1713244923791704];
         data.locs = [-0.9324695142031521
                -0.6612093864662645
                -0.2386191860831969
                0.2386191860831969
                0.6612093864662645
                0.9324695142031521];
    case 7
        wgts = [0.417959184
                0.381830051
                0.381830051
                0.279705391
                0.279705391
                0.129484966
                0.129484966 ];
            data.locs = [0
                    0.405845151
                    -0.405845151
                    -0.741531186
                    0.741531186
                    -0.949107912
                    0.949107912
                    ];
    case 8
        wgts = [0.362683783
                0.362683783
                0.313706646
                0.313706646
                0.222381034
                0.222381034
                0.101228536
                0.101228536
                ];
        data.locs = [-0.183434642
                0.183434642
                -0.52553241
                0.52553241
                -0.796666477
                0.796666477
                -0.960289856
                0.960289856
                ];
    case 9
        wgts = [0.330239355
                0.180648161
                0.180648161
                0.081274388
                0.081274388
                0.312347077
                0.312347077
                0.260610696
                0.260610696
                ];
        data.locs = [0
                -0.836031107
                0.836031107
                -0.96816024
                0.96816024
                -0.324253423
                0.324253423
                -0.613371433
                0.613371433
                ];
    case 10
        wgts = [0.295524225
                0.295524225
                0.269266719
                0.269266719
                0.219086363
                0.219086363
                0.149451349
                0.149451349
                0.066671344
                0.066671344
                ];
        data.locs = [-0.148874339
                0.148874339
                -0.433395394
                0.433395394
                -0.679409568
                0.679409568
                -0.865063367
                0.865063367
                -0.973906529
                0.973906529
                ];
    case 11       
        wgts = [0.272925087
                0.262804545
                0.262804545
                0.233193765
                0.233193765
                0.186290211
                0.186290211
                0.125580369
                0.125580369
                0.055668567
                0.055668567
                ];
        data.locs = [0
                -0.269543156
                0.269543156
                -0.519096129
                0.519096129
                -0.730152006
                0.730152006
                -0.8870626
                0.8870626
                -0.978228658
                0.978228658
                ];
    case 20
        mat = [1	0.1527533871307258	-0.0765265211334973
            2	0.1527533871307258	0.0765265211334973
            3	0.1491729864726037	-0.2277858511416451
            4	0.1491729864726037	0.2277858511416451
            5	0.1420961093183820	-0.3737060887154195
            6	0.1420961093183820	0.3737060887154195
            7	0.1316886384491766	-0.5108670019508271
            8	0.1316886384491766	0.5108670019508271
            9	0.1181945319615184	-0.6360536807265150
            10	0.1181945319615184	0.6360536807265150
            11	0.1019301198172404	-0.7463319064601508
            12	0.1019301198172404	0.7463319064601508
            13	0.0832767415767048	-0.8391169718222188
            14	0.0832767415767048	0.8391169718222188
            15	0.0626720483341091	-0.9122344282513259
            16	0.0626720483341091	0.9122344282513259
            17	0.0406014298003869	-0.9639719272779138
            18	0.0406014298003869	0.9639719272779138
            19	0.0176140071391521	-0.9931285991850949
            20	0.0176140071391521	0.9931285991850949];
        wgts = mat(:,2);
        data.locs = mat(:,3);
    case 64       
        wgts = [            
            0.048690957
            0.048690957
            0.048575467
            0.048575467
            0.048344762
            0.048344762
            0.047999389
            0.047999389
            0.047540166
            0.047540166
            0.046968183
            0.046968183
            0.046284797
            0.046284797
            0.045491628
            0.045491628
            0.044590558
            0.044590558
            0.043583725
            0.043583725
            0.042473515
            0.042473515
            0.041262563
            0.041262563
            0.039953741
            0.039953741
            0.038550153
            0.038550153
            0.037055129
            0.037055129
            0.035472213
            0.035472213
            0.033805162
            0.033805162
            0.032057928
            0.032057928
            0.030234657
            0.030234657
            0.028339673
            0.028339673
            0.02637747
            0.02637747
            0.024352703
            0.024352703
            0.022270174
            0.022270174
            0.020134823
            0.020134823
            0.017951716
            0.017951716
            0.01572603
            0.01572603
            0.013463048
            0.013463048
            0.011168139
            0.011168139
            0.00884676
            0.00884676
            0.006504458
            0.006504458
            0.004147033
            0.004147033
            0.001783281
            0.001783281];
        data.locs = [-0.024350293
                0.024350293
                -0.072993122
                0.072993122
                -0.121462819
                0.121462819
                -0.16964442
                0.16964442
                -0.217423644
                0.217423644
                -0.264687162
                0.264687162
                -0.311322872
                0.311322872
                -0.357220158
                0.357220158
                -0.402270158
                0.402270158
                -0.446366017
                0.446366017
                -0.489403146
                0.489403146
                -0.531279464
                0.531279464
                -0.571895646
                0.571895646
                -0.611155355
                0.611155355
                -0.648965471
                0.648965471
                -0.685236313
                0.685236313
                -0.71988185
                0.71988185
                -0.752819907
                0.752819907
                -0.783972359
                0.783972359
                -0.813265315
                0.813265315
                -0.840629296
                0.840629296
                -0.865999398
                0.865999398
                -0.889315446
                0.889315446
                -0.910522137
                0.910522137
                -0.929569172
                0.929569172
                -0.946411375
                0.946411375
                -0.9610088
                0.9610088
                -0.973326828
                0.973326828
                -0.983336254
                0.983336254
                -0.991013371
                0.991013371
                -0.996340117
                0.996340117
                -0.999305042
                0.999305042
                ];
    otherwise
            disp('Error')
end
data.P = length(data.locs);
data.wgts1D = wgts; 

data.n1quad = data.P;
data.n2quad = 1;
data.n3quad = 1;
if data.space_dims >= 2
    data.n2quad = data.P;
    if data.space_dims >= 3
        data.n3quad = data.P;
    end
end

data.n1limiter = data.P+2;
data.n2limiter = 1;
data.n3limiter = 1;
if data.space_dims >= 2
    data.n2limiter = data.P+2;
    if data.space_dims >= 3
        data.n3limiter = data.P+2;
    end
end

switch data.space_dims
    case 1
        for v1quad = 1:data.n1quad
            data.wgts_space(v1quad,1) = wgts(v1quad);
            for tquad = 1:data.n1quad
                data.wgts_spacetime(tquad,v1quad) = wgts(tquad)*wgts(v1quad);
            end
        end
    case 2
        for v1quad = 1:data.n1quad
        for v2quad = 1:data.n2quad
            data.wgts_space(v1quad,v2quad) = wgts(v1quad)*wgts(v2quad);
            for tquad = 1:data.n1quad
                data.wgts_spacetime(tquad,v1quad,v2quad) = wgts(tquad)*wgts(v1quad)*wgts(v2quad);
            end
        end
        end
    case 3
        for v1quad = 1:data.n1quad
        for v2quad = 1:data.n2quad
        for v3quad = 1:data.n3quad
            data.wgts_space(v1quad,v2quad,v3quad) = wgts(v1quad)*wgts(v2quad)*wgts(v3quad);
            for tquad = 1:data.n1quad
                data.wgts_spacetime(tquad,v1quad,v2quad,v3quad) = wgts(tquad)*wgts(v1quad)*wgts(v2quad)*wgts(v3quad);
            end
        end
        end
        end
end

end


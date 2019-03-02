function [data] = initialize_method_implicit(data)
% written by Pierson Guthrey

switch data.predictorbasis
    case 'Q'
        switch data.space_dims
        %M (order) , r+1 (parameter) , D (dimensions)
        %--------------------------------------------------
        %              1D CASE
        %--------------------------------------------------
            case 1
            %LI 1D case
            cfl_restrictions(1,1,1)= 1;
            cfl_restrictions(2,1,1)= .3; %??????
            cfl_restrictions(3,1,1)= .1;
            cfl_restrictions(4,1,1)= .104;
            cfl_restrictions(5,1,1)= .06;
            cfl_restrictions(6,1,1)= .04;%??????
            cfl_restrictions(8,1,1)= .03;
            cfl_restrictions(10,1,1)= .01;

            %RI1 1D
            cfl_restrictions(1,2,1)= .9;
            cfl_restrictions(2,2,1)= .9;
            cfl_restrictions(3,2,1)= .9;
            cfl_restrictions(4,2,1)= .9;
            cfl_restrictions(5,2,1)= .9;
            cfl_restrictions(6,2,1)= .9;
            cfl_restrictions(8,2,1)= .9;
            cfl_restrictions(10,2,1)= .9;
            cfl_restrictions(20,2,1)= .9;

            %RI2 1D
            cfl_restrictions(1,3,1)= 1.2;
            cfl_restrictions(2,3,1)= 1.2;
            cfl_restrictions(4,3,1)= 1.2;%1.2;
            cfl_restrictions(6,3,1)= 1.2;
            cfl_restrictions(8,3,1)= 1.2;
            cfl_restrictions(10,3,1)= 1.2;

            %RI3 1D
            cfl_restrictions(4,4,1)= 1.2;

        %--------------------------------------------------
        %              2D CASE
        %--------------------------------------------------    
            case 2

            %LI 2D case
            cfl_restrictions(1,1,2)= .8;
            cfl_restrictions(2,1,2)= .1;%.2;
            cfl_restrictions(3,1,2)= .1;
            cfl_restrictions(4,1,2)= .05;%.06;
            cfl_restrictions(6,1,2)= .03;%.06; % .06 < cfl < .7
            cfl_restrictions(8,1,2)= .02;
            cfl_restrictions(10,1,2)= .01;

            %RI1 2D
            cfl_restrictions(1,2,2)= 1;
            cfl_restrictions(2,2,2)= .80;%4; %cfl > .6 ?????
            cfl_restrictions(3,2,2)= .75;%4; %cfl > .6 ?????
            cfl_restrictions(4,2,2)= .75;%.75;
            cfl_restrictions(6,2,2)= .75;%.75;%.5; %  < cfl < .5
            cfl_restrictions(8,2,2)= .7;%.75;%.5; %  < cfl < .5
            cfl_restrictions(10,2,2)= .75;%.5; %  < cfl < .5

            %RI2 2D
            cfl_restrictions(1,3,2)= 1;
            cfl_restrictions(4,3,2)= 1;
            cfl_restrictions(6,3,2)= 1;%.7; %  < cfl < .9
            cfl_restrictions(10,3,2)= 1;%.7; %  < cfl < .9

            %RI3 2D
            cfl_restrictions(4,4,2)= 1;
            cfl_restrictions(6,4,2)= 1;%.7; %  < cfl < .9

        %--------------------------------------------------
        %              3D CASE
        %--------------------------------------------------
            case 3
            %LI 3D case
            cfl_restrictions(1,1,3)= 1;
            cfl_restrictions(2,1,3)= .1;
            cfl_restrictions(4,1,3)= .03;
            cfl_restrictions(6,1,3)= .025;
            cfl_restrictions(8,1,3)= .02;
            cfl_restrictions(10,1,3)= .01;

            %RI1 3D
            cfl_restrictions(1,2,3)= 1;%.5; %  < cfl < .5
            cfl_restrictions(2,2,3)= .8;%.5; %  < cfl < .5
            cfl_restrictions(4,2,3)= .6;%75;%.5; %  < cfl < .5
            cfl_restrictions(6,2,3)= .6;%75;%.5; %  < cfl < .5
            cfl_restrictions(8,2,3)= .6;%75;%.5; %  < cfl < .5
            cfl_restrictions(10,2,3)= .6;%.5; %  < cfl < .5

            %RI2 3D
            cfl_restrictions(1,3,3)= 1;%.7; %  < cfl < .9
            cfl_restrictions(4,3,3)= 1;%.7; %  < cfl < .9
        end
    case 'P' 
        switch data.space_dims
        %M (order) , r+1 (parameter) , D (dimensions)
        %--------------------------------------------------
        %              1D CASE
        %--------------------------------------------------
            case 1

            %LI 1D case
            cfl_restrictions(1,1,1)= 1;
            cfl_restrictions(2,1,1)= .3; %??????
            cfl_restrictions(3,1,1)= .1;
            cfl_restrictions(4,1,1)= .104;
            cfl_restrictions(5,1,1)= .06;
            cfl_restrictions(6,1,1)= .04;%??????
            cfl_restrictions(8,1,1)= .03;
            cfl_restrictions(10,1,1)= .01;

            %RI1 1D
            cfl_restrictions(1,2,1)= .9;
            cfl_restrictions(2,2,1)= .9;
            cfl_restrictions(3,2,1)= .9;
            cfl_restrictions(4,2,1)= .9;
            cfl_restrictions(5,2,1)= .9;
            cfl_restrictions(6,2,1)= .9;
            cfl_restrictions(8,2,1)= .9;
            cfl_restrictions(10,2,1)= 1;
            cfl_restrictions(20,2,1)= 1;

            %RI2 1D
            cfl_restrictions(1,3,1)= 1.2;
            cfl_restrictions(2,3,1)= 1.2;
            cfl_restrictions(4,3,1)= 1.2;%1.2;
            cfl_restrictions(6,3,1)= 1.2;
            cfl_restrictions(8,3,1)= 1.2;
            cfl_restrictions(10,3,1)= 1.2;

            %RI3 1D
            cfl_restrictions(4,4,1)= 1.2;

        %--------------------------------------------------
        %              2D CASE
        %--------------------------------------------------    
            case 2

            %LI 2D case
            cfl_restrictions(1,1,2)= .8;
            cfl_restrictions(2,1,2)= .1;%.2;
            cfl_restrictions(3,1,2)= .1;
            cfl_restrictions(4,1,2)= .05;%.06;
            cfl_restrictions(6,1,2)= .03;%.06; % .06 < cfl < .7
            cfl_restrictions(8,1,2)= .02;
            cfl_restrictions(10,1,2)= .01;

            %RI1 2D
            cfl_restrictions(1,2,2)= 1;
            cfl_restrictions(2,2,2)= .80;%4; %cfl > .6 ?????
            cfl_restrictions(3,2,2)= .7;%4; %cfl > .6 ?????
            cfl_restrictions(4,2,2)= .7;%.75;
            cfl_restrictions(6,2,2)= .7;%.75;%.5; %  < cfl < .5
            cfl_restrictions(8,2,2)= .7;%.5; %  < cfl < .5
            cfl_restrictions(10,2,2)= .7;%.5; %  < cfl < .5

            %RI2 2D
            cfl_restrictions(1,3,2)= 1;
            cfl_restrictions(4,3,2)= 1;
            cfl_restrictions(6,3,2)= 1;%.7; %  < cfl < .9
            cfl_restrictions(10,3,2)= 1;%.7; %  < cfl < .9

            %RI3 2D
            cfl_restrictions(4,4,2)= 1;
            cfl_restrictions(6,4,2)= 1;%.7; %  < cfl < .9

        %--------------------------------------------------
        %              3D CASE
        %--------------------------------------------------
            case 3
            %LI 3D case
            cfl_restrictions(1,1,3)= 1;
            cfl_restrictions(2,1,3)= .1;
            cfl_restrictions(4,1,3)= .03;
            cfl_restrictions(6,1,3)= .025;
            cfl_restrictions(8,1,3)= .02;
            cfl_restrictions(10,1,3)= .01;

            %RI1 3D
            cfl_restrictions(1,2,3)= 1;%.5; %  < cfl < .5
            cfl_restrictions(2,2,3)= .8;%.5; %  < cfl < .5
            cfl_restrictions(4,2,3)= .6;%75;%.5; %  < cfl < .5
            cfl_restrictions(6,2,3)= .6;%75;%.5; %  < cfl < .5
            cfl_restrictions(8,2,3)= .6;%75;%.5; %  < cfl < .5
            cfl_restrictions(10,2,3)= .6;%.5; %  < cfl < .5

            %RI2 3D
            cfl_restrictions(1,3,3)= 1;%.7; %  < cfl < .9
            cfl_restrictions(4,3,3)= 1;%.7; %  < cfl < .9        
        end
end
data.cfl = cfl_restrictions(data.M,data.r_param+1,data.space_dims);

end       
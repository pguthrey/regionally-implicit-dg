function [fcoeffs_new] = opsplit_Phase(deltat,fcoeffs,E1coeffs,A1coeffs,A2coeffs,data)

data.Tfinal = deltat;

%h1 = waitbar(0,'Advecting in data.Phase','OuterPosition', [800 100 300 75]);
h1 = waitbar(0,'Phase Space Operator...','OuterPosition', [100 400 300 75]);
normphi = 8;

for xquad = 1:data.P
for v1quad = 1:data.P
for v2quad = 1:data.P
    zetatrans = data.zeta1D2V(1,:,xquad,v1quad,v2quad)';
    data.PHI = data.vectphi(1,:,v1quad,v2quad);
    wgt = data.wgts3D(xquad,v1quad,v2quad);
    ProjM(:,:,xquad) = zetatrans*data.PHI*wgt/8;
end
end
end

fcoeffs_new = zeros(size(fcoeffs));

for ix = 1:data.Nx
    proj = zeros(data.theta1D2V,1);
    for xquad = 1:data.P
        xii = data.locs(xquad); 
        x = data.xcenters(ix)+xii*data.deltax/2;
        A1 = DGeval_1D(A1coeffs,x,data);
        A2 = DGeval_1D(A2coeffs,x,data);
        data.appdata.x = x;
        data.appdata.E1 = DGeval_1D(E1coeffs,x,data);
        data.appdata.E2 = A1 + A2;
        data.appdata.B3 = A1 - A2;    
        
        data.appdata.fcoeffs = fcoeffs;
        DG_initialconditions = projection_DGL2_Proj(@(point) problem_IC(point,data),data);
        
        tnow = 0;
        nstep = 2;
        data.Time = 0;
        check = 1;
        deltatnew = deltat;
        DGsolution_old = DG_initialconditions;

        while check 
            timeleft = data.Tfinal - tnow;
            if deltat > timeleft
                %Reduce timestep
                if data.verbose
                    disp('final step, trying time left as deltat')
                end
                deltat = timeleft;
                data.nuv1 = deltat/data.deltav1;
                data.nuv2 = deltat/data.deltav2;
                data.nuv3 = deltat/data.deltav3;
            end    
            if data.verbose
                disp(['Attempting dt = ' num2str(deltat) ' at t = ' num2str(tnow)])
            end
  
            [DGprediction,predmaxF,predmaxG,predmaxH] = predictor_main_implicit(DGsolution_old,data);
            DGprediction_limited = DGprediction;   
            [DGcorrected,corrmaxF,corrmaxG,corrmaxH] = corrector_main_noghosts(DGprediction_limited,DGsolution_old,data); 
            data.deltat = deltat;
            DGsolution_new = limiter_phi(DGcorrected,DGsolution_old,data);
            maxspeedx = max([predmaxF corrmaxF]);
            maxspeedy = max([predmaxG corrmaxG]);
            maxspeedz = max([predmaxH corrmaxH]);
            dtspeedF = data.deltav1/maxspeedx*data.cfl;
            dtspeedG = data.deltav2/maxspeedy*data.cfl;
            dtspeedH = data.deltav3/maxspeedz*data.cfl;
            deltatnew = min([dtspeedF dtspeedG dtspeedH]);
            if deltat > deltatnew
                %Reject timestep
                if data.verbose
                    disp('reject step, trying smaller dt')
                end
                deltat = .95*deltatnew;% deltatnew*.95;
                data.nuv1 = deltat/data.deltav1;
                data.nuv2 = deltat/data.deltav2;    
                data.nuv3 = deltat/data.deltav3;    
            elseif (deltat == deltatnew)
                %Accept Timestep, no changes to dt
                if data.verbose
                    disp('accept, keep same dt')
                end
                tnow = tnow + deltat;
                data.Time(nstep) = tnow ; 
                if data.makegif_conserved
                    [data] = plotting_makegif(DGsolution_new,tnow,data,nstep);
                elseif data.plotwhilerunning 
                    problem_plotting_whilerunning(DGsolution_new,tnow,data)
                end
                nstep = nstep + 1;
                DGsolution_old = DGsolution_new;
                done = tnow/data.Tfinal;
                waitbar(done, h1);
            else 
                %Accept timestep, increase dt
                if data.verbose
                    disp('accept, change dt')
                end
                tnow = tnow + deltat;
                deltat = .95*deltatnew;
                data.nuv1 = deltat/data.deltav1;
                data.nuv2 = deltat/data.deltav2;
                data.nuv3 = deltat/data.deltav3;
                data.Time(nstep) = tnow ;
                if data.makegif_conserved
                    [data] = plotting_makegif(DGsolution_new,tnow,data,nstep);
                elseif data.plotwhilerunning 
                    problem_plotting_whilerunning(DGsolution_new,tnow,data)
                end
                nstep = nstep + 1;
                DGsolution_old = DGsolution_new;
                %done = tnow/data.Tfinal;
                %waitbar(done, h1);
            end
            %disp('------------------------------------')
            check = (tnow < data.Tfinal);
        end
        delete(h1)
        DGsolution_end = DGsolution_new;
        for iv1 = 1:data.Nv1
        for iv2 = 1:data.Nv2            
            coeffs(:,1) = DGsolution_end(:,iv1,iv2);
            fcoeffs_new(:,ix,iv1,iv2) = fcoeffs_new(:,ix,iv1,iv2) + ProjM(:,:,xquad)*coeffs;
        end
        end
    end
    %waitbar(ix/Nx,h1)
end
%keyboard



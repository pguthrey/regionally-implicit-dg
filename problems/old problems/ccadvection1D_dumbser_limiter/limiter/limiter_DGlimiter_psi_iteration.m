function [ DGcoeffs_limited ] = limiter_DGlimiter_psi(DGcoeffs,DGpast,data)

DGcoeffs_limited = DGcoeffs;
debug_plot_DGpsi(DGcoeffs_limited,data)
disp('pre limiting')
pause

xmesh = 1:data.Nv1;
DGcoeffs_limited = DGcoeffs;
check = 1;
iters = 0;
[ BETA ] = limiter_DGlimiter_detector(DGcoeffs_limited,DGpast,data);
while check

    iters = iters + 1
    
    for iv1  = xmesh(BETA==1)
        Averages = data.V*DGcoeffs_limited(:,iv1)
        
        
        DGcoeffs_limited(:,iv1) = data.limiter_filter*;
    end
    
    [ BETA ] = limiter_DGlimiter_detector(DGcoeffs_limited,DGpast,data);




    check = any(BETA == 1) && (iters < 1000);
end

    debug_plot_DGpsi(DGcoeffs_limited,data)
    disp('post limiting')
    pause
keyboard
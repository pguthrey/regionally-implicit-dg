function [searchdir] = predictor_RI1DG_solver(Jacobian,residual,thetaT)%#codegen
% written by Pierson Guthrey

    inds = @(l,r) (1+ (l-1)*thetaT ):( thetaT*r) ;
    ind  = @(i) (1:thetaT) + (i-1)*thetaT;
    
    %ROW #1 ------------------------------------------------------------
    this_row = 1;
    this_cols = [2,4];
    
    block_this_row = ind(this_row);
    block_this_cols = inds(this_cols(1),this_cols(2));
    
    
    blockLHS = Jacobian(block_this_row,block_this_row);      
    Jacobian(block_this_row,block_this_cols) = blockLHS\Jacobian(block_this_row,block_this_cols);
    residual(block_this_row) = blockLHS\residual(block_this_row);
    
    target_row = 2;    
    block_target_row = ind(target_row);    
    tempblock = Jacobian(block_target_row,block_this_row);
    Jacobian(block_target_row,block_this_cols) = Jacobian(block_target_row,block_this_cols) - tempblock*Jacobian(block_this_row,block_this_cols);
    residual(block_target_row) = residual(block_target_row) - tempblock*residual(block_this_row);
    
    target_row = 4;    
    block_target_row = ind(target_row);    
    tempblock = Jacobian(block_target_row,block_this_row);
    Jacobian(block_target_row,block_this_cols) = Jacobian(block_target_row,block_this_cols) - tempblock*Jacobian(block_this_row,block_this_cols);
    residual(block_target_row) = residual(block_target_row) - tempblock*residual(block_this_row);
    
    %ROW #2 ------------------------------------------------------------
    this_row = 2;
    this_cols = [3,5];
    
    block_this_row = ind(this_row);
    block_this_cols = inds(this_cols(1),this_cols(2));
    
    
    blockLHS = Jacobian(block_this_row,block_this_row);      
    Jacobian(block_this_row,block_this_cols) = blockLHS\Jacobian(block_this_row,block_this_cols);
    residual(block_this_row) = blockLHS\residual(block_this_row);
    
    target_row = 3;    
    block_target_row = ind(target_row);    
    tempblock = Jacobian(block_target_row,block_this_row);
    Jacobian(block_target_row,block_this_cols) = Jacobian(block_target_row,block_this_cols) - tempblock*Jacobian(block_this_row,block_this_cols);
    residual(block_target_row) = residual(block_target_row) - tempblock*residual(block_this_row);
    
    target_row = 4;    
    block_target_row = ind(target_row);    
    tempblock = Jacobian(block_target_row,block_this_row);
    Jacobian(block_target_row,block_this_cols) = Jacobian(block_target_row,block_this_cols) - tempblock*Jacobian(block_this_row,block_this_cols);
    residual(block_target_row) = residual(block_target_row) - tempblock*residual(block_this_row);
    
    target_row = 5;    
    block_target_row = ind(target_row);    
    tempblock = Jacobian(block_target_row,block_this_row);
    Jacobian(block_target_row,block_this_cols) = Jacobian(block_target_row,block_this_cols) - tempblock*Jacobian(block_this_row,block_this_cols);
    residual(block_target_row) = residual(block_target_row) - tempblock*residual(block_this_row);
    
    %ROW #3 ------------------------------------------------------------
    this_row = 3;
    this_cols = [4,6];
    
    block_this_row = ind(this_row);
    block_this_cols = inds(this_cols(1),this_cols(2));
    
    
    blockLHS = Jacobian(block_this_row,block_this_row);      
    Jacobian(block_this_row,block_this_cols) = blockLHS\Jacobian(block_this_row,block_this_cols);
    residual(block_this_row) = blockLHS\residual(block_this_row);
        
    target_row = 4;    
    block_target_row = ind(target_row);    
    tempblock = Jacobian(block_target_row,block_this_row);
    Jacobian(block_target_row,block_this_cols) = Jacobian(block_target_row,block_this_cols) - tempblock*Jacobian(block_this_row,block_this_cols);
    residual(block_target_row) = residual(block_target_row) - tempblock*residual(block_this_row);
    
    target_row = 5;    
    block_target_row = ind(target_row);    
    tempblock = Jacobian(block_target_row,block_this_row);
    Jacobian(block_target_row,block_this_cols) = Jacobian(block_target_row,block_this_cols) - tempblock*Jacobian(block_this_row,block_this_cols);
    residual(block_target_row) = residual(block_target_row) - tempblock*residual(block_this_row);
        
    target_row = 6;    
    block_target_row = ind(target_row);    
    tempblock = Jacobian(block_target_row,block_this_row);
    Jacobian(block_target_row,block_this_cols) = Jacobian(block_target_row,block_this_cols) - tempblock*Jacobian(block_this_row,block_this_cols);
    residual(block_target_row) = residual(block_target_row) - tempblock*residual(block_this_row);

    %ROW #4 ------------------------------------------------------------
    this_row = 4;
    this_cols = [5,7];
    
    block_this_row = ind(this_row);
    block_this_cols = inds(this_cols(1),this_cols(2));
    
    
    blockLHS = Jacobian(block_this_row,block_this_row);      
    Jacobian(block_this_row,block_this_cols) = blockLHS\Jacobian(block_this_row,block_this_cols);
    residual(block_this_row) = blockLHS\residual(block_this_row);
        
    target_row = 5;    
    block_target_row = ind(target_row);    
    tempblock = Jacobian(block_target_row,block_this_row);
    Jacobian(block_target_row,block_this_cols) = Jacobian(block_target_row,block_this_cols) - tempblock*Jacobian(block_this_row,block_this_cols);
    residual(block_target_row) = residual(block_target_row) - tempblock*residual(block_this_row);    
    
    target_row = 6;    
    block_target_row = ind(target_row);    
    tempblock = Jacobian(block_target_row,block_this_row);
    Jacobian(block_target_row,block_this_cols) = Jacobian(block_target_row,block_this_cols) - tempblock*Jacobian(block_this_row,block_this_cols);
    residual(block_target_row) = residual(block_target_row) - tempblock*residual(block_this_row);    
    
    target_row = 7;    
    block_target_row = ind(target_row);    
    tempblock = Jacobian(block_target_row,block_this_row);
    Jacobian(block_target_row,block_this_cols) = Jacobian(block_target_row,block_this_cols) - tempblock*Jacobian(block_this_row,block_this_cols);
    residual(block_target_row) = residual(block_target_row) - tempblock*residual(block_this_row);    
    
    %ROW #9 ------------------------------------------------------------
    this_row = 9;
    this_cols = [6,8];
    
    block_this_row = ind(this_row);
    block_this_cols = inds(this_cols(1),this_cols(2));
    
    
    blockLHS = Jacobian(block_this_row,block_this_row);      
    Jacobian(block_this_row,block_this_cols) = blockLHS\Jacobian(block_this_row,block_this_cols);
    residual(block_this_row) = blockLHS\residual(block_this_row);
        
    target_row = 8;    
    block_target_row = ind(target_row);    
    tempblock = Jacobian(block_target_row,block_this_row);
    Jacobian(block_target_row,block_this_cols) = Jacobian(block_target_row,block_this_cols) - tempblock*Jacobian(block_this_row,block_this_cols);
    residual(block_target_row) = residual(block_target_row) - tempblock*residual(block_this_row);    
    
    target_row = 6;    
    block_target_row = ind(target_row);    
    tempblock = Jacobian(block_target_row,block_this_row);
    Jacobian(block_target_row,block_this_cols) = Jacobian(block_target_row,block_this_cols) - tempblock*Jacobian(block_this_row,block_this_cols);
    residual(block_target_row) = residual(block_target_row) - tempblock*residual(block_this_row);    
    
    %ROW #8 ------------------------------------------------------------
    this_row = 8;
    this_cols = [5,7];
    
    block_this_row = ind(this_row);
    block_this_cols = inds(this_cols(1),this_cols(2));
    
    
    blockLHS = Jacobian(block_this_row,block_this_row);      
    Jacobian(block_this_row,block_this_cols) = blockLHS\Jacobian(block_this_row,block_this_cols);
    residual(block_this_row) = blockLHS\residual(block_this_row);
        
    target_row = 7;    
    block_target_row = ind(target_row);    
    tempblock = Jacobian(block_target_row,block_this_row);
    Jacobian(block_target_row,block_this_cols) = Jacobian(block_target_row,block_this_cols) - tempblock*Jacobian(block_this_row,block_this_cols);
    residual(block_target_row) = residual(block_target_row) - tempblock*residual(block_this_row);    
    
    target_row = 6;    
    block_target_row = ind(target_row);    
    tempblock = Jacobian(block_target_row,block_this_row);
    Jacobian(block_target_row,block_this_cols) = Jacobian(block_target_row,block_this_cols) - tempblock*Jacobian(block_this_row,block_this_cols);
    residual(block_target_row) = residual(block_target_row) - tempblock*residual(block_this_row);    
    
    target_row = 5;    
    block_target_row = ind(target_row);    
    tempblock = Jacobian(block_target_row,block_this_row);
    Jacobian(block_target_row,block_this_cols) = Jacobian(block_target_row,block_this_cols) - tempblock*Jacobian(block_this_row,block_this_cols);
    residual(block_target_row) = residual(block_target_row) - tempblock*residual(block_this_row);    
    
    %ROW #7 ------------------------------------------------------------
    this_row = 7;
    this_cols = [5,6];
    
    block_this_row = ind(this_row);
    block_this_cols = inds(this_cols(1),this_cols(2));
    
    
    blockLHS = Jacobian(block_this_row,block_this_row);      
    Jacobian(block_this_row,block_this_cols) = blockLHS\Jacobian(block_this_row,block_this_cols);
    residual(block_this_row) = blockLHS\residual(block_this_row);
        
    target_row = 6;    
    block_target_row = ind(target_row);    
    tempblock = Jacobian(block_target_row,block_this_row);
    Jacobian(block_target_row,block_this_cols) = Jacobian(block_target_row,block_this_cols) - tempblock*Jacobian(block_this_row,block_this_cols);
    residual(block_target_row) = residual(block_target_row) - tempblock*residual(block_this_row);    
    
    target_row = 5;    
    block_target_row = ind(target_row);    
    tempblock = Jacobian(block_target_row,block_this_row);
    Jacobian(block_target_row,block_this_cols) = Jacobian(block_target_row,block_this_cols) - tempblock*Jacobian(block_this_row,block_this_cols);
    residual(block_target_row) = residual(block_target_row) - tempblock*residual(block_this_row);    
    
    %ROW #6 ------------------------------------------------------------
    this_row = 6;
    this_cols = [5,5];
    
    block_this_row = ind(this_row);
    block_this_cols = inds(this_cols(1),this_cols(2));    
    
    blockLHS = Jacobian(block_this_row,block_this_row);      
    Jacobian(block_this_row,block_this_cols) = blockLHS\Jacobian(block_this_row,block_this_cols);
    residual(block_this_row) = blockLHS\residual(block_this_row);
        
    target_row = 5;    
    block_target_row = ind(target_row);    
    tempblock = Jacobian(block_target_row,block_this_row);
    Jacobian(block_target_row,block_this_cols) = Jacobian(block_target_row,block_this_cols) - tempblock*Jacobian(block_this_row,block_this_cols);
    residual(block_target_row) = residual(block_target_row) - tempblock*residual(block_this_row);   
    
    % Finalize
    this_row = 5;    
    block_this_row = ind(this_row);    
    blockLHS = Jacobian(block_this_row,block_this_row);      
    residual(block_this_row) = blockLHS\residual(block_this_row);
    
    searchdir = residual(block_this_row);

end


function [Jacobian,residual] = problem_region_Jacobian_residual_M2_P(q,q_past,nuv1) 
Upsi_1_3 = q(3,1,1);
Upsi_2_3 = q(3,2,1);
Upsi_3_3 = q(3,3,1);
Upsi_1_2 = q(2,1,1);
Upsi_2_2 = q(2,2,1);
Upsi_3_2 = q(2,3,1);
Upsi_1_1 = q(1,1,1);
Upsi_2_1 = q(1,2,1);
Upsi_3_1 = q(1,3,1);
Upsi_past_1_3 = q_past(3,1,1);
Upsi_past_2_3 = q_past(3,2,1);
Upsi_past_3_3 = q_past(3,3,1);
Upsi_past_1_2 = q_past(2,1,1);
Upsi_past_2_2 = q_past(2,2,1);
Upsi_past_3_2 = q_past(2,3,1);
Upsi_past_1_1 = q_past(1,1,1);
Upsi_past_2_1 = q_past(1,2,1);
Upsi_past_3_1 = q_past(1,3,1);
Uphi1_east_2 = Upsi_1_3;
Uphi2_east_2 = Upsi_2_3;
Uphi3_east_2 = Upsi_3_3;
Uphi1_east_1 = Upsi_1_1 + 1.73205080756887729*Upsi_1_2;
Uphi2_east_1 = Upsi_2_1 + 1.73205080756887729*Upsi_2_2;
Uphi3_east_1 = Upsi_3_1 + 1.73205080756887729*Upsi_3_2;
Uphi1_west_2 = Upsi_1_3;
Uphi2_west_2 = Upsi_2_3;
Uphi3_west_2 = Upsi_3_3;
Uphi1_west_1 = Upsi_1_1 - 1.73205080756887729*Upsi_1_2;
Uphi2_west_1 = Upsi_2_1 - 1.73205080756887729*Upsi_2_2;
Uphi3_west_1 = Upsi_3_1 - 1.73205080756887729*Upsi_3_2;
lambda_sample_12_1 = max(abs(Uphi1_east_1 - 1.0*Uphi1_east_2),abs(Uphi2_west_1 - 1.0*Uphi2_west_2));
lambda_sample_12_2 = max(abs(Uphi1_east_1 + Uphi1_east_2),abs(Uphi2_west_1 + Uphi2_west_2));
lambda_sample_23_1 = max(abs(Uphi2_east_1 - 1.0*Uphi2_east_2),abs(Uphi3_west_1 - 1.0*Uphi3_west_2));
lambda_sample_23_2 = max(abs(Uphi2_east_1 + Uphi2_east_2),abs(Uphi3_west_1 + Uphi3_west_2));
lambda_12_1 = 0.5*lambda_sample_12_1_1 + 0.5*lambda_sample_12_2_1;
lambda_12_2 = 0.5*lambda_sample_12_2_1 - 0.5*lambda_sample_12_1_1;
lambda_23_1 = 0.5*lambda_sample_23_1_1 + 0.5*lambda_sample_23_2_1;
lambda_23_2 = 0.5*lambda_sample_23_2_1 - 0.5*lambda_sample_23_1_1;
residual = zeros(9,1);
residual(9,1) = 6.0*Upsi_3_3 - 3.46410161513775459*Upsi_3_1 + 3.46410161513775459*Upsi_past_3_1 + 6.0*Upsi_past_3_3 + nuv1*(2.0*Uphi3_east_1*Uphi3_east_2 - 1.0*Uphi2_east_1*Uphi2_east_2 - 1.0*Uphi3_west_1*Uphi3_west_2 - 1.0*Uphi2_east_1*lambda_23_2 - 1.0*Uphi2_east_2*lambda_23_1 + Uphi3_west_1*lambda_23_2 + Uphi3_west_2*lambda_23_1);
residual(8,1) = 2.0*Upsi_3_2 - 2.0*Upsi_past_3_2 + nuv1*(1.73205080756887729*Uphi2_east_1*lambda_23_1 + 1.73205080756887729*Uphi2_east_2*lambda_23_2 - 1.73205080756887729*Uphi3_west_1*lambda_23_1 - 1.73205080756887729*Uphi3_west_2*lambda_23_2 + 0.866025403784438647*Uphi2_east_1^2 + 0.866025403784438647*Uphi2_east_2^2 + 1.73205080756887729*Uphi3_east_1^2 + 1.73205080756887729*Uphi3_east_2^2 + 0.866025403784438647*Uphi3_west_1^2 + 0.866025403784438647*Uphi3_west_2^2 - 3.46410161513775459*Upsi_3_1^2 - 3.46410161513775459*Upsi_3_2^2 - 3.46410161513775459*Upsi_3_3^2);
residual(7,1) = 2.0*Upsi_3_1 + 3.46410161513775459*Upsi_3_3 - 2.0*Upsi_past_3_1 - 3.46410161513775459*Upsi_past_3_3 + nuv1*(Uphi3_west_1*lambda_23_1 - 1.0*Uphi2_east_2*lambda_23_2 - 1.0*Uphi2_east_1*lambda_23_1 + Uphi3_west_2*lambda_23_2 - 0.5*Uphi2_east_1^2 - 0.5*Uphi2_east_2^2 + Uphi3_east_1^2 + Uphi3_east_2^2 - 0.5*Uphi3_west_1^2 - 0.5*Uphi3_west_2^2);
residual(6,1) = 6.0*Upsi_2_3 - 3.46410161513775459*Upsi_2_1 + 3.46410161513775459*Upsi_past_2_1 + 6.0*Upsi_past_2_3 + nuv1*(Uphi2_east_1*Uphi2_east_2 - 1.0*Uphi1_east_1*Uphi1_east_2 - 1.0*Uphi2_west_1*Uphi2_west_2 + Uphi3_west_1*Uphi3_west_2 - 1.0*Uphi1_east_1*lambda_12_2 - 1.0*Uphi1_east_2*lambda_12_1 + Uphi2_west_1*lambda_12_2 + Uphi2_west_2*lambda_12_1 + Uphi2_east_1*lambda_23_2 + Uphi2_east_2*lambda_23_1 - 1.0*Uphi3_west_1*lambda_23_2 - 1.0*Uphi3_west_2*lambda_23_1);
residual(5,1) = 2.0*Upsi_2_2 - 2.0*Upsi_past_2_2 + nuv1*(1.73205080756887729*Uphi1_east_1*lambda_12_1 + 1.73205080756887729*Uphi1_east_2*lambda_12_2 - 1.73205080756887729*Uphi2_west_1*lambda_12_1 - 1.73205080756887729*Uphi2_west_2*lambda_12_2 + 1.73205080756887729*Uphi2_east_1*lambda_23_1 + 1.73205080756887729*Uphi2_east_2*lambda_23_2 - 1.73205080756887729*Uphi3_west_1*lambda_23_1 - 1.73205080756887729*Uphi3_west_2*lambda_23_2 + 0.866025403784438647*Uphi1_east_1^2 + 0.866025403784438647*Uphi1_east_2^2 + 0.866025403784438647*Uphi2_east_1^2 + 0.866025403784438647*Uphi2_east_2^2 + 0.866025403784438647*Uphi2_west_1^2 + 0.866025403784438647*Uphi2_west_2^2 + 0.866025403784438647*Uphi3_west_1^2 + 0.866025403784438647*Uphi3_west_2^2 - 3.46410161513775459*Upsi_2_1^2 - 3.46410161513775459*Upsi_2_2^2 - 3.46410161513775459*Upsi_2_3^2);
residual(4,1) = 2.0*Upsi_2_1 + 3.46410161513775459*Upsi_2_3 - 2.0*Upsi_past_2_1 - 3.46410161513775459*Upsi_past_2_3 + nuv1*(Uphi2_west_1*lambda_12_1 - 1.0*Uphi1_east_2*lambda_12_2 - 1.0*Uphi1_east_1*lambda_12_1 + Uphi2_west_2*lambda_12_2 + Uphi2_east_1*lambda_23_1 + Uphi2_east_2*lambda_23_2 - 1.0*Uphi3_west_1*lambda_23_1 - 1.0*Uphi3_west_2*lambda_23_2 - 0.5*Uphi1_east_1^2 - 0.5*Uphi1_east_2^2 + 0.5*Uphi2_east_1^2 + 0.5*Uphi2_east_2^2 - 0.5*Uphi2_west_1^2 - 0.5*Uphi2_west_2^2 + 0.5*Uphi3_west_1^2 + 0.5*Uphi3_west_2^2);
residual(3,1) = 6.0*Upsi_1_3 - 3.46410161513775459*Upsi_1_1 + 3.46410161513775459*Upsi_past_1_1 + 6.0*Upsi_past_1_3 + nuv1*(Uphi1_east_1*Uphi1_east_2 - 2.0*Uphi1_west_1*Uphi1_west_2 + Uphi2_west_1*Uphi2_west_2 + Uphi1_east_1*lambda_12_2 + Uphi1_east_2*lambda_12_1 - 1.0*Uphi2_west_1*lambda_12_2 - 1.0*Uphi2_west_2*lambda_12_1);
residual(2,1) = 2.0*Upsi_1_2 - 2.0*Upsi_past_1_2 + nuv1*(1.73205080756887729*Uphi1_east_1*lambda_12_1 + 1.73205080756887729*Uphi1_east_2*lambda_12_2 - 1.73205080756887729*Uphi2_west_1*lambda_12_1 - 1.73205080756887729*Uphi2_west_2*lambda_12_2 + 0.866025403784438647*Uphi1_east_1^2 + 0.866025403784438647*Uphi1_east_2^2 + 1.73205080756887729*Uphi1_west_1^2 + 1.73205080756887729*Uphi1_west_2^2 + 0.866025403784438647*Uphi2_west_1^2 + 0.866025403784438647*Uphi2_west_2^2 - 3.46410161513775459*Upsi_1_1^2 - 3.46410161513775459*Upsi_1_2^2 - 3.46410161513775459*Upsi_1_3^2);
residual(1,1) = 2.0*Upsi_1_1 + 3.46410161513775459*Upsi_1_3 - 2.0*Upsi_past_1_1 - 3.46410161513775459*Upsi_past_1_3 + nuv1*(Uphi1_east_1*lambda_12_1 + Uphi1_east_2*lambda_12_2 - 1.0*Uphi2_west_1*lambda_12_1 - 1.0*Uphi2_west_2*lambda_12_2 + 0.5*Uphi1_east_1^2 + 0.5*Uphi1_east_2^2 - 1.0*Uphi1_west_1^2 - 1.0*Uphi1_west_2^2 + 0.5*Uphi2_west_1^2 + 0.5*Uphi2_west_2^2);
Jacobian = zeros(9,9);
Jacobian(9,9) = 6.0 + nuv1*(2.0*Uphi3_east_1 - 1.0*Uphi3_west_1 + lambda_23_1);
Jacobian(9,8) =  + nuv1*(3.46410161513775459*Uphi3_east_2 + 1.73205080756887729*Uphi3_west_2 - 1.73205080756887729*lambda_23_2);
Jacobian(9,7) = -3.46410161513775459 + nuv1*(2.0*Uphi3_east_2 - 1.0*Uphi3_west_2 + lambda_23_2);
Jacobian(9,6) =  + nuv1*(- 1.0*Uphi2_east_1 - 1.0*lambda_23_1);
Jacobian(9,5) =  + nuv1*(- 1.73205080756887729*Uphi2_east_2 - 1.73205080756887729*lambda_23_2);
Jacobian(9,4) =  + nuv1*(- 1.0*Uphi2_east_2 - 1.0*lambda_23_2);
Jacobian(8,9) =  + nuv1*(3.46410161513775459*Uphi3_east_2 + 1.73205080756887729*Uphi3_west_2 - 6.92820323027550917*Upsi_3_3 - 1.73205080756887729*lambda_23_2);
Jacobian(8,8) = 2.0 + nuv1*(6.0*Uphi3_east_1 - 3.0*Uphi3_west_1 - 6.92820323027550917*Upsi_3_2 + 3.0*lambda_23_1);
Jacobian(8,7) =  + nuv1*(3.46410161513775459*Uphi3_east_1 + 1.73205080756887729*Uphi3_west_1 - 6.92820323027550917*Upsi_3_1 - 1.73205080756887729*lambda_23_1);
Jacobian(8,6) =  + nuv1*(1.73205080756887729*Uphi2_east_2 + 1.73205080756887729*lambda_23_2);
Jacobian(8,5) =  + nuv1*(3.0*Uphi2_east_1 + 3.0*lambda_23_1);
Jacobian(8,4) =  + nuv1*(1.73205080756887729*Uphi2_east_1 + 1.73205080756887729*lambda_23_1);
Jacobian(7,9) = 3.46410161513775459 + nuv1*(2.0*Uphi3_east_2 - 1.0*Uphi3_west_2 + lambda_23_2);
Jacobian(7,8) =  + nuv1*(3.46410161513775459*Uphi3_east_1 + 1.73205080756887729*Uphi3_west_1 - 1.73205080756887729*lambda_23_1);
Jacobian(7,7) = 2.0 + nuv1*(2.0*Uphi3_east_1 - 1.0*Uphi3_west_1 + lambda_23_1);
Jacobian(7,6) =  + nuv1*(- 1.0*Uphi2_east_2 - 1.0*lambda_23_2);
Jacobian(7,5) =  + nuv1*(- 1.73205080756887729*Uphi2_east_1 - 1.73205080756887729*lambda_23_1);
Jacobian(7,4) =  + nuv1*(- 1.0*Uphi2_east_1 - 1.0*lambda_23_1);
Jacobian(6,9) =  + nuv1*(Uphi3_west_1 - 1.0*lambda_23_1);
Jacobian(6,8) =  + nuv1*(1.73205080756887729*lambda_23_2 - 1.73205080756887729*Uphi3_west_2);
Jacobian(6,7) =  + nuv1*(Uphi3_west_2 - 1.0*lambda_23_2);
Jacobian(6,6) = 6.0 + nuv1*(Uphi2_east_1 - 1.0*Uphi2_west_1 + lambda_12_1 + lambda_23_1);
Jacobian(6,5) =  + nuv1*(1.73205080756887729*Uphi2_east_2 + 1.73205080756887729*Uphi2_west_2 - 1.73205080756887729*lambda_12_2 + 1.73205080756887729*lambda_23_2);
Jacobian(6,4) = -3.46410161513775459 + nuv1*(Uphi2_east_2 - 1.0*Uphi2_west_2 + lambda_12_2 + lambda_23_2);
Jacobian(6,3) =  + nuv1*(- 1.0*Uphi1_east_1 - 1.0*lambda_12_1);
Jacobian(6,2) =  + nuv1*(- 1.73205080756887729*Uphi1_east_2 - 1.73205080756887729*lambda_12_2);
Jacobian(6,1) =  + nuv1*(- 1.0*Uphi1_east_2 - 1.0*lambda_12_2);
Jacobian(5,9) =  + nuv1*(1.73205080756887729*Uphi3_west_2 - 1.73205080756887729*lambda_23_2);
Jacobian(5,8) =  + nuv1*(3.0*lambda_23_1 - 3.0*Uphi3_west_1);
Jacobian(5,7) =  + nuv1*(1.73205080756887729*Uphi3_west_1 - 1.73205080756887729*lambda_23_1);
Jacobian(5,6) =  + nuv1*(1.73205080756887729*Uphi2_east_2 + 1.73205080756887729*Uphi2_west_2 - 6.92820323027550917*Upsi_2_3 - 1.73205080756887729*lambda_12_2 + 1.73205080756887729*lambda_23_2);
Jacobian(5,5) = 2.0 + nuv1*(3.0*Uphi2_east_1 - 3.0*Uphi2_west_1 - 6.92820323027550917*Upsi_2_2 + 3.0*lambda_12_1 + 3.0*lambda_23_1);
Jacobian(5,4) =  + nuv1*(1.73205080756887729*Uphi2_east_1 + 1.73205080756887729*Uphi2_west_1 - 6.92820323027550917*Upsi_2_1 - 1.73205080756887729*lambda_12_1 + 1.73205080756887729*lambda_23_1);
Jacobian(5,3) =  + nuv1*(1.73205080756887729*Uphi1_east_2 + 1.73205080756887729*lambda_12_2);
Jacobian(5,2) =  + nuv1*(3.0*Uphi1_east_1 + 3.0*lambda_12_1);
Jacobian(5,1) =  + nuv1*(1.73205080756887729*Uphi1_east_1 + 1.73205080756887729*lambda_12_1);
Jacobian(4,9) =  + nuv1*(Uphi3_west_2 - 1.0*lambda_23_2);
Jacobian(4,8) =  + nuv1*(1.73205080756887729*lambda_23_1 - 1.73205080756887729*Uphi3_west_1);
Jacobian(4,7) =  + nuv1*(Uphi3_west_1 - 1.0*lambda_23_1);
Jacobian(4,6) = 3.46410161513775459 + nuv1*(Uphi2_east_2 - 1.0*Uphi2_west_2 + lambda_12_2 + lambda_23_2);
Jacobian(4,5) =  + nuv1*(1.73205080756887729*Uphi2_east_1 + 1.73205080756887729*Uphi2_west_1 - 1.73205080756887729*lambda_12_1 + 1.73205080756887729*lambda_23_1);
Jacobian(4,4) = 2.0 + nuv1*(Uphi2_east_1 - 1.0*Uphi2_west_1 + lambda_12_1 + lambda_23_1);
Jacobian(4,3) =  + nuv1*(- 1.0*Uphi1_east_2 - 1.0*lambda_12_2);
Jacobian(4,2) =  + nuv1*(- 1.73205080756887729*Uphi1_east_1 - 1.73205080756887729*lambda_12_1);
Jacobian(4,1) =  + nuv1*(- 1.0*Uphi1_east_1 - 1.0*lambda_12_1);
Jacobian(3,6) =  + nuv1*(Uphi2_west_1 - 1.0*lambda_12_1);
Jacobian(3,5) =  + nuv1*(1.73205080756887729*lambda_12_2 - 1.73205080756887729*Uphi2_west_2);
Jacobian(3,4) =  + nuv1*(Uphi2_west_2 - 1.0*lambda_12_2);
Jacobian(3,3) = 6.0 + nuv1*(Uphi1_east_1 - 2.0*Uphi1_west_1 + lambda_12_1);
Jacobian(3,2) =  + nuv1*(1.73205080756887729*Uphi1_east_2 + 3.46410161513775459*Uphi1_west_2 + 1.73205080756887729*lambda_12_2);
Jacobian(3,1) = -3.46410161513775459 + nuv1*(Uphi1_east_2 - 2.0*Uphi1_west_2 + lambda_12_2);
Jacobian(2,6) =  + nuv1*(1.73205080756887729*Uphi2_west_2 - 1.73205080756887729*lambda_12_2);
Jacobian(2,5) =  + nuv1*(3.0*lambda_12_1 - 3.0*Uphi2_west_1);
Jacobian(2,4) =  + nuv1*(1.73205080756887729*Uphi2_west_1 - 1.73205080756887729*lambda_12_1);
Jacobian(2,3) =  + nuv1*(1.73205080756887729*Uphi1_east_2 + 3.46410161513775459*Uphi1_west_2 - 6.92820323027550917*Upsi_1_3 + 1.73205080756887729*lambda_12_2);
Jacobian(2,2) = 2.0 + nuv1*(3.0*Uphi1_east_1 - 6.0*Uphi1_west_1 - 6.92820323027550917*Upsi_1_2 + 3.0*lambda_12_1);
Jacobian(2,1) =  + nuv1*(1.73205080756887729*Uphi1_east_1 + 3.46410161513775459*Uphi1_west_1 - 6.92820323027550917*Upsi_1_1 + 1.73205080756887729*lambda_12_1);
Jacobian(1,6) =  + nuv1*(Uphi2_west_2 - 1.0*lambda_12_2);
Jacobian(1,5) =  + nuv1*(1.73205080756887729*lambda_12_1 - 1.73205080756887729*Uphi2_west_1);
Jacobian(1,4) =  + nuv1*(Uphi2_west_1 - 1.0*lambda_12_1);
Jacobian(1,3) = 3.46410161513775459 + nuv1*(Uphi1_east_2 - 2.0*Uphi1_west_2 + lambda_12_2);
Jacobian(1,2) =  + nuv1*(1.73205080756887729*Uphi1_east_1 + 3.46410161513775459*Uphi1_west_1 + 1.73205080756887729*lambda_12_1);
Jacobian(1,1) = 2.0 + nuv1*(Uphi1_east_1 - 2.0*Uphi1_west_1 + lambda_12_1);

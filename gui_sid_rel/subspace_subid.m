function [A B C D X0 SSM] = subspace_subid(u,y,time_points, n)
% n is the order of model(A is n x n matrix)
% u is input in time domain [x1,x2,...], x_i has the same length of
% time_points

% y is output in time domain
% time_points is time vector  (must be linear)

% A,B, C and D are system matrices

addpath ../gui_sid/vanoverschee/
order_number = n;
%the following 2 lines are useless for sid_det
n_port = size(u,1); 
n_blk_row = 2*ceil(n/n_port);
n_blk_row = order_number+1;
Ts = time_points(2)-time_points(1);  %  Ts is a scalar value equal to the sampling interval of your experiment.

[A,B,C,D, K, R] = subid(y,u,n_blk_row);

%sys_sid_det = sid_det(y',u',order_number,'alg1');
%A = sys_sid_det.A;
%B = sys_sid_det.B;
%C = sys_sid_det.C;
%D = sys_sid_det.D;
% K=zeros(size(A,1),size(C,1)); 
X0=zeros(size(A(1,:)));
%X0(1)=1;



SSM = idss(A, B, C, D, K, X0',Ts);


function [A B C D X0 SSM] = subspace_pem(u,y,time_points, n)
% n is the order of model(A is n x n matrix)
% u is input in time domain [x1,x2,...], x_i has the same length of
% time_points

% y is output in time domain
% time_points is time vector  (must be linear)

% A,B, C and D are system matrices


order_number = n;
%the following 2 lines are useless for sid_det
n_port = size(u,1); 
n_blk_row = 2*ceil(n/n_port);

Ts = time_points(2)-time_points(1);  %  Ts is a scalar value equal to the sampling interval of your experiment.

%[A,B,C,D] = subid(y,u,n_blk_row);
%K=zeros(size(A,1),size(C,1)); 
%x0=zeros(size(A(1,:)));

%sys_sid_det = sid_det(y',u',order_number,'alg1');
%A = sys_sid_det.A;
%B = sys_sid_det.B;
%C = sys_sid_det.C;
%D = sys_sid_det.D;
%K=zeros(size(A,1),size(C,1)); 
%X0=zeros(size(A(1,:)));


%SSM = idss(A, B, C, D, K, X0',Ts);
data = iddata(y,u,Ts);   % set up input data format
SSM = pem(data,order_number,'focus','Simulation','cov','None','InitialState','estimate');    % you can choose other model u like, for here, it looks better than pem(stable) 
%SSM = n4sid(data,order_number,'focus','Simulation','cov','None','InitialState','estimate');
%SSM = n4sid(data,order_number,'focus','Simulation','N4Weight','CVA','cov','None','InitialState','estimate','n4h',[ 32 11 11 ]);
%SSM = n4sid(data,order_number,'focus','simulation','N4Weight','Auto','cov','None','InitialState','estimate','N4Horizon','auto');
%SSM = n4sid(data,order_number,'focus','Simulation','N4Weight','CVA','cov','None','InitialState','estimate','n4h',[30 10 10],'nk',zeros(1,size(u,2)));
%SSM = n4sid(data,order_number,'focus','Simulation','N4Weight','CVA','cov','None','InitialState','estimate','n4h',[27 10 10]);
%SSM = n4sid(data,order_number,'focus','Simulation','N4Weight','CVA','cov','None','InitialState','estimate','N4Horizon',[40 10 10],'nk',[ones(1,size(u,2))]);
%SSM = n4sid(data,order_number,'focus','Simulation','N4Weight','CVA','cov','None','InitialState','estimate','n4h',[40 10 10]);


A = SSM.A;
B = SSM.B;
C = SSM.C;
D = SSM.D;
X0 = SSM.X0;

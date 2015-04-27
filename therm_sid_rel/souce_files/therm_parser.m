function therm_parser
%parses excel data and extracts power inputs and temperature outputs for
%each core.  Outputs data in .mat files

H = xlsread('QuadCore_Al_UCR.xls','Training_Al_UCR');
H2 = xlsread('QuadCore_Al_UCR.xls','Training_Al_UCR (2)');
input = xlsread('QuadCore_Al_UCR.xls','Training_Al_UCR_powers');
input2 = xlsread('QuadCore_Al_UCR.xls','Training_Al_UCR_powers (2)');

time = H(:,1);
die_cen = H(:,2);
s0 = H(:,3);
s1 = H(:,4);
s2 = H(:,5);
s3 = H(:,6);
s4 = H(:,7);
HH = H(:,8);
II = H(:,9);
JJ = H(:,10);
KK = H(:,11);
LL = H(:,12);
MM = H(:,13);
NN = H(:,14);
OO = H(:,15);
PP = H(:,16);
QQ = H(:,17);
RR = H(:,18);

time_v = H2(:,1);
die_cen_v = H2(:,2);
s0_v = H2(:,3);
s1_v = H2(:,4);
s2_v = H2(:,5);
s3_v = H2(:,6);
s4_v = H2(:,7);
HH_v = H(:,8);
II_v = H(:,9);
JJ_v = H(:,10);
KK_v = H(:,11);
LL_v = H(:,12);
MM_v = H(:,13);
NN_v = H(:,14);
OO_v = H(:,15);
PP_v = H(:,16);
QQ_v = H(:,17);
RR_v = H(:,18);

time_in = input(:,1);
s0_in = input(:,2);
s1_in = input(:,3);
s2_in = input(:,4);
s3_in = input(:,5);
s4_in = input(:,6);

time_in_v = input2(:,1);
s0_in_v = input2(:,2);
s1_in_v = input2(:,3);
s2_in_v = input2(:,4);
s3_in_v = input2(:,5);
s4_in_v = input2(:,6);

die_cen_u = zeros(length(s4_in),1);
s0_u = zeros(length(s4_in),1);
s1_u = zeros(length(s4_in),1);
s2_u = zeros(length(s4_in),1);
s3_u = zeros(length(s4_in),1);
s4_u = zeros(length(s4_in),1);

die_cen_u_v = zeros(length(s4_in),1);
s0_u_v = zeros(length(s4_in),1);
s1_u_v = zeros(length(s4_in),1);
s2_u_v = zeros(length(s4_in),1);
s3_u_v = zeros(length(s4_in),1);
s4_u_v = zeros(length(s4_in),1);

for k = 1:length(time_in)-1
    
    for j = 1:length(time)
        
        if( (time_in(k) <= time(j)) && (time(j) < time_in(k+1)) )
            s0_u(j) = s0_in(k);
            s1_u(j) = s1_in(k);
            s2_u(j) = s2_in(k);
            s3_u(j) = s3_in(k);
            s4_u(j) = s4_in(k);
            
        end
        
        if( (time_in_v(k) <= time(j)) && (time(j) < time_in_v(k+1)) )
            s0_u_v(j) = s0_in_v(k);
            s1_u_v(j) = s1_in_v(k);
            s2_u_v(j) = s2_in_v(k);
            s3_u_v(j) = s3_in_v(k);
            s4_u_v(j) = s4_in_v(k);            
        end
        
    end
end

s0_u(length(time)) = s0_in(length(time_in));
s1_u(length(time)) = s1_in(length(time_in));
s2_u(length(time)) = s2_in(length(time_in));
s3_u(length(time)) = s3_in(length(time_in));
s4_u(length(time)) = s4_in(length(time_in));

s0_u_v(length(time)) = s0_in_v(length(time_in));
s1_u_v(length(time)) = s1_in_v(length(time_in));
s2_u_v(length(time)) = s2_in_v(length(time_in));
s3_u_v(length(time)) = s3_in_v(length(time_in));
s4_u_v(length(time)) = s4_in_v(length(time_in));

savefile = 'therm_data.mat';
save(savefile,'die_cen','s0','s1','s2','s3','s4','s0_u','s1_u','s2_u','s3_u','s4_u','time','HH','II','JJ','KK','LL','MM','NN','OO','PP','QQ','RR');
savefile = 'therm_data_v.mat';
save(savefile,'s0_in_v','s1_in_v','s2_in_v','s3_in_v','s4_in_v','die_cen_v','s0_v','s1_v','s2_v','s3_v','s4_v','s0_u_v','s1_u_v','s2_u_v','s3_u_v','s4_u_v','time_v','HH_v','II_v','JJ_v','KK_v','LL_v','MM_v','NN_v','OO_v','PP_v','QQ_v','RR_v');
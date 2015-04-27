%function error = cross_valid(uu, yy, time, order, partitions)

clear all;
clc;
partitions = 10;
order = 25;

% load therm_data.mat;
% load therm_data_v.mat;
load therm_datanew_par.mat;
load therm_datanew_v.mat;

load data_orig.mat;

%time = times;
Ts = time(2) - time(1);
% %yy = [ die_cen s0 s1 s2 s3 s4 ];
  yy = [ s0 s1 s2 s3 s4 ];
  uu = [ s0_u s1_u s2_u s3_u s4_u ];
  
  %data_orig = iddata(yy,uu,Ts);
% 
yy_v = [ s0_v s1_v s2_v s3_v s4_v ];
uu_v = [ s0_u_v s1_u_v s2_u_v s3_u_v s4_u_v ];

%yy = [ yy' ex_yy' ]';
%uu = [ uu' ex_uu' ]';

% yy = yy0;
% uu = uu0;

% data_var2 = iddata(yy1,uu1,Ts);
% data_var1 = iddata(yy2,uu2,Ts);

%time = times;

data_var = iddata(yy_v, uu_v, Ts);
data = iddata(yy, uu, Ts);

p_index = round(length(uu)/partitions);

first_index = 1;
last_index = p_index;
len = length(uu);

fit_er = [];
SSMM = cell(1,partitions);

for j = 1:partitions
    
   train_uu = [];
   train_yy = [];
      
   valid_uu = uu(first_index:last_index,:);
   valid_yy = yy(first_index:last_index,:);
   
   if( first_index > 1 )
    train_uu = uu(1:first_index,:);
    train_yy = yy(1:first_index,:);
   end
   
   train_uu = [ train_uu' uu(last_index:len,:)' ]';
   train_yy = [ train_yy' yy(last_index:len,:)' ]'; 
   
   data_train = iddata(train_yy, train_uu, Ts);
   data_valid = iddata(valid_yy, valid_uu, Ts);
   [A B C D X0 SSM] = subspace(train_uu,train_yy,time,order);
   
   if( j == 9 )
       uu_hold = train_uu;
       yy_hold = train_yy;
   end
   
   SSMM{1,j} = SSM;
   
   [yh,fit,x0] = compare(data_var,SSM);
   %figure; compare(data_valid,SSM);
   %[yh,fitt,x0] = compare(data_train,SSM);
   %[yh,fit,x0] = compare(data,SSM);
   fit_er = [ fit_er mean(fit) ];
   
%    figure
%    compare(data_train,SSM);
  
   first_index = last_index;
   check_index = last_index;
   last_index = last_index + p_index;
   
   if( last_index > length(uu) || last_index == length(uu) )
       last_index = length(uu);
   end
   
end

[ A B C D X0 SSM ] = subspace(uu,yy,time,order);
%[ A B C D X0 SSM2 ] = subspace(uu,yy,time,order);

savefile = 'cross_val.mat';
save(savefile,'SSMM');

%load data_orig.mat;
fitter = [];
fir = []
for j = 1:partitions
    [yh,fit,x0] = compare(data,SSMM{1,j});
    [yh,fitt,x0] = compare(data,SSM);    
    %[yh,fitt,x0] = compare(data_orig,SSMM{1,j});
    fitter = [ fitter [mean(fit); mean(fitt)] ];
    fir(j) = (mean(fit)+fit_er(j))/2;
end

fitty = [ fitter; fit_er ];

max_ind = [];
%fit_err = fit_er;
fit_err = fir;
for j = 1:5
    [C,max_ind(j)] = max(fit_err);
    fit_err(max_ind(j)) = 0;
end

yy_new = [];
uu_new = [];

for j = 1:length(yy)
    if( j >= max_ind(1)*50-50 && j <= max_ind(1)*50 )
    elseif( j >= max_ind(3)*50-50 && j <= max_ind(3)*50 )
    else
        yy_new = [ yy_new; yy(j,:) ];
        uu_new = [ uu_new; uu(j,:) ];
    end
end
% 
% [A B C D X0 SSX] = subspace(uu_new,yy_new,time,order);
% [yh,fit,x0] = compare(data,SSX);
% mean(fit)
   
% figure;
% compare(data_orig, SSM);
% figure;
% compare(data_orig, SSMM{1,max_ind(1)});
% figure;
% compare(data_var, SSM);
% figure;
% compare(data_var, SSMM{1,max_ind(1)}); 
% figure;
% compare(data_orig, SSX);

%data_orig = data_orig1;

% data_ap = iddata(yy0(1:800,:),uu0(1:800,:),times(2)-times(1));
% data_ex = data_ap;

%[y_ssm1,fit,x0] = compare(data_orig, SSM);
%[y_ssmm1,fit,x0] = compare(data_orig, SSMM{1,max_ind(1)});
[y_ssm1,fit,x0] = compare(data_orig, SSM);
[y_ssmm1,fit,x0] = compare(data_orig, SSMM{1,max_ind(1)});
% [y_ssm2,fit,x0] = compare(data_var, SSM);
% [y_ssmm2,fit,x0] = compare(data_var, SSMM{1,max_ind(1)}); 

% [y_ssm1,fit,x0] = compare(data_orig, SSM);
% [y_ssm2,fit,x0] = compare(data_orig, SSM2);

%[data_er,data_er2] = cross_valid_plot(yy_ex,y_ssm1,y_ssmm1,time_ex);
[data_er,data_er2] = cross_valid_plot(yy0,y_ssm1,y_ssmm1,times);
%[data_er,data_er2] = cross_valid_plot(data_orig.OutputData(:,:),y_ssm1,y_ssmm1,times);

% figure;
% plot(times,yy0(:,3),'k','linewidth',2); hold on;
% plot(times,y_ssm1{1,1}.OutputData(:,3),':r','linewidth',2); hold on;
% plot(times,y_ssm2{1,1}.OutputData(:,3),'--g','linewidth',2); hold on;
% legend('verif','data set 1','data set 2');
% xlabel('time (sec)');
% ylabel('temp (celcius)');
% title('Temperature Response Core 2');



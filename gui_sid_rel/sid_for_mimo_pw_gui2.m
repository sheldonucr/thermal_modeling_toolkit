clc;
clear all;
addpath /data/zliu/
%load siso_temp_pw.mat;
load prbs_mimo_2.mat;
%load prbs_siso6_2.mat;
app_start = 1;
order = 4;%3;%15;%5;
app_end = 2*3400;%5100;%3200;
%app_len1 = 3400;%3200;
app_len = app_end;%4800;
app_tot = 20400;
app_offset = 1700;%1700; %1600;
offset = 293.15;
ver_start=1;
ver_end = size(y_temp_vec,1);
if(mod(app_tot,app_len))
    stop;
else
NumPw = app_tot/(app_len-app_offset);
end

Ts = t(2)-t(1);

beg = 1;

y_app1 = y_temp_vec(app_start:app_end,:)'- offset;
u_app1 = p_vec(:,app_start:app_end);
[ A1 B1 C1 D1 X01 SSM1 ] = subspace(u_app1', y_app1', t(app_start:app_end)', order);
xh1(:,1)=X01;
for i=1:app_len
  xh1(:,i+1) = A1*xh1(:,i) + B1*p_vec(:,i);
  yh(:,i) = C1*xh1(:,i) + D1*p_vec(:,i);
end
app_start = app_start + app_len - app_offset ;
app_end = app_start + app_len-1;
iter1=0;

% 
% T_pre = eye(size(A1));
% xh_pre= xh1(:,end-app_offset:end-1);
% for k=2:11
% y_app = y_temp_vec(app_start:app_end,:)'- offset;
% u_app = p_vec(:,app_start:app_end);
% [ A B C D X0 SSM ] = subspace(u_app', y_app', t(app_start:app_end)', order);
% xh(:,1)=X0;
% for i=1:app_len
%   xh(:,i+1) = A*xh(:,i) + B*p_vec(:,i+app_start-1);
%   yh(:,i+app_start-1) = C*xh(:,i) + D*p_vec(:,i+app_start-1);
% end
% T = xh(:,1:app_offset)/xh_pre*T_pre;
% T_inverse = inv(T);
% 
% At = T_inverse*A*T;
% Bt = T_inverse*B;
% Ct = C*T;
% Dt = D;
% T_pre = T;
% xh_pre=xh(:,end-app_offset:end-1);
% app_start = app_start + app_len - app_offset;
% app_end = app_start + app_len-1;
% fid=fopen('ss_models.m','w');
% Nk = k;
% fprintf(fid, 'SSM%d=SSM;\n',Nk);
% fprintf(fid, 'A%dt=At; B%dt=Bt; C%dt=Ct; D%dt=Dt;\n',Nk,Nk,Nk,Nk);
% fprintf(fid, 'disp(''generate A%dt, B%dt, C%dt, D%dt'')',Nk,Nk,Nk,Nk);
% fclose(fid);
% Nk
% run ss_models;
%delete ss_models.*
%end

fid=fopen('ss_models.m','w');
T_pre = eye(size(A1));
xh_pre= xh1(:,end-app_offset:end-1);
for k=2:NumPw-1
y_app = y_temp_vec(app_start:app_end,:)'- offset;
u_app = p_vec(:,app_start:app_end);
[ A(:,:,k) B(:,:,k) C(:,:,k) D(:,:,k) X0(:,k) SSM ] = subspace(u_app', y_app', t(app_start:app_end)', order);
xh(:,1)=X0(:,k);
for i=1:app_len
  xh(:,i+1) = reshape(A(:,:,k),size(A1))*xh(:,i) + reshape(B(:,:,k),size(B1))*p_vec(:,i+app_start-1);
  yh(:,i+app_start-1) = reshape(C(:,:,k),size(C1))*xh(:,i) + reshape(D(:,:,k),size(D1))*p_vec(:,i+app_start-1);
end

T = xh(:,1:app_offset)/xh_pre*T_pre;
T_inverse = inv(T);

At(:,:,k) = T_inverse*reshape(A(:,:,k), size(A1))*T;
Bt(:,:,k) = T_inverse*reshape(B(:,:,k),size(B1));
Ct(:,:,k) = reshape(C(:,:,k),size(C1))*T;
Dt(:,:,k) = reshape(D(:,:,k),size(D1));
T_pre = T;
xh_pre=xh(:,end-app_offset:end-1);
if(k<NumPw-1)
app_start = app_start + app_len - app_offset;
app_end = app_start + app_len-1;
end

fprintf(fid, 'A%dt=reshape(At(:,:,%d),size(A1)); B%dt=reshape(Bt(:,:,%d),size(B1)); C%dt=reshape(Ct(:,:,%d),size(C1)); D%dt=reshape(Dt(:,:,%d),size(D1));iter%d=0;\n',k,k,k,k,k,k,k,k,k);
end

%check the remaining data  and decide whether or not to build another model
if(app_end ~= app_tot)
    k=k+1;
    app_start = app_start + app_len - app_offset;
    y_app = y_temp_vec(app_start:app_tot,:)'- offset;
    u_app = p_vec(:,app_start:app_tot);
    [ A(:,:,k) B(:,:,k) C(:,:,k) D(:,:,k) X0(:,k) SSM ] = subspace(u_app', y_app', t(app_start:app_tot)', order);
    xh(:,1)=X0(:,k);

    for i=1:app_tot-app_end+app_offset
    xh(:,i+1) = reshape(A(:,:,k),size(A1))*xh(:,i) + reshape(B(:,:,k),size(B1))*p_vec(:,i+app_start-1);
    yh(:,i+app_start-1) = reshape(C(:,:,k),size(C1))*xh(:,i) + reshape(D(:,:,k),size(D1))*p_vec(:,i+app_start-1);
    end
    T = xh(:,1:app_offset)/xh_pre*T_pre;
    T_inverse = inv(T);

    At(:,:,k) = T_inverse*reshape(A(:,:,k), size(A1))*T;
    Bt(:,:,k) = T_inverse*reshape(B(:,:,k),size(B1));
    Ct(:,:,k) = reshape(C(:,:,k),size(C1))*T;
    Dt(:,:,k) = reshape(D(:,:,k),size(D1));
   
    fprintf(fid, 'A%dt=reshape(At(:,:,%d),size(A1)); B%dt=reshape(Bt(:,:,%d),size(B1)); C%dt=reshape(Ct(:,:,%d),size(C1)); D%dt=reshape(Dt(:,:,%d),size(D1));iter%d=0;\n',k,k,k,k,k,k,k,k,k);

end
%end
    
app_end=app_tot;
fclose(fid);
run ss_models;


y_ver = y_temp_vec(app_start:ver_end,:)' - offset;
u_ver = p_vec(:,app_start:ver_end)';

%data_ver = iddata(y_ver,u_ver,Ts);

figure; plot(t(1:ver_end)',y_temp_vec(1:ver_end,1)'-273.15 ,'g','linewidth',2); hold on;

 xh=zeros(order,ver_end);
 xh(:,1) = X01;
 yh=zeros(25,ver_end);
 yh(:,1) = C1*xh(:,1) + D1*p_vec(:,1);
 %iter1=0; iter2=0; iter3=0; iter4=0; iter5=0; iter6=0; iter7=0; iter8=0; iter9=0; iter10=0; iter11=0; 
 %setup the simulation data partition
 TempTransitPoint = 32;
 TempStep = 20;
 
 fid2=fopen('sim_model.m','w');
 fprintf(fid2, 'if(yh_m<%f)\n',TempTransitPoint);
 fprintf(fid2,'xh(:,i+1) = A1*xh(:,i) + B1*p_vec(:,i);\n');
 fprintf(fid2,'yh(:,i) = C1*xh(:,i) + D1*p_vec(:,i);\n');
 fprintf(fid2,'if(i>app_end) iter%d=iter%d+1; end\n',1,1);
 for m=2:k %k is the number of the model
     if(m<k)
      fprintf(fid2, 'elseif(yh_m>%f && yh_m<%f)\n',TempTransitPoint+10*(m-2),TempTransitPoint+TempStep*(m-1));
     else
      fprintf(fid2, 'elseif(yh_m>%f)\n',TempTransitPoint+TempStep*(m-2));  
     end     
      fprintf(fid2,'xh(:,i+1) = A%dt*xh(:,i) + B%dt*p_vec(:,i);\n',m,m);
      fprintf(fid2,'yh(:,i) = C%dt*xh(:,i) + D%dt*p_vec(:,i);\n',m,m);
      fprintf(fid2,'if(i>app_end) iter%d=iter%d+1; end\n',m,m);
 end
 fprintf(fid2,'end');
 fclose(fid2);
 tic
 for i=1:ver_end
 if(i==1)
   yh_m = sum(yh(:,1))/25+20;
 else
     yh_m = y_temp_vec(i)-offset+20; %sum(yh(:,i-1))/25+20; %y_temp_vec(i)-offset+20; %yh(i-1);%y_temp_vec(i);
 end

  run sim_model
%   if(yh_m < 32)
%    xh(:,i+1) = A1*xh(:,i) + B1*p_vec(:,i);
%    yh(:,i) = C1*xh(:,i) + D1*p_vec(:,i);
%    if(i>app_end) iter1=iter1+1; end
%   elseif ( (yh_m> 32) && (yh_m <36))
%    xh(:,i+1) = A2t*xh(:,i) + B2t*p_vec(:,i);
%    yh(:,i) = C2t*xh(:,i) + D2t*p_vec(:,i);
%    if(i>app_end) iter2=iter2+1; end
%   elseif ( (yh_m> 36) && (yh_m <41))
%    xh(:,i+1) = A3t*xh(:,i) + B3t*p_vec(:,i);
%    yh(:,i) = C3t*xh(:,i) + D3t*p_vec(:,i);
%    if(i>app_end) iter3=iter3+1; end
%   elseif ((yh_m> 41) && (yh_m < 46))
%    xh(:,i+1) = A4t*xh(:,i) + B4t*p_vec(:,i);
%    yh(:,i) = C4t*xh(:,i) + D4t*p_vec(:,i);
%    if(i>app_end) iter4=iter4+1; end
%   elseif ((yh_m> 46) && (yh_m < 51.5))
%    xh(:,i+1) = A5t*xh(:,i) + B5t*p_vec(:,i);
%    yh(:,i) = C5t*xh(:,i) + D5t*p_vec(:,i);
%    if(i>app_end) iter5=iter5+1; end
%   elseif ((yh_m> 51.5) && (yh_m < 56.5))
%    xh(:,i+1) = A6t*xh(:,i) + B6t*p_vec(:,i);
%    yh(:,i) = C6t*xh(:,i) + D6t*p_vec(:,i); 
%     if(i>app_end) iter6=iter6+1; end
%   elseif ((yh_m>56.5) && (yh_m<61.5))
%    xh(:,i+1) = A7t*xh(:,i) + B7t*p_vec(:,i);
%    yh(:,i) = C7t*xh(:,i) + D7t*p_vec(:,i);  
%     if(i>app_end) iter7=iter7+1; end  
%   elseif ((yh_m>61.5) && (yh_m<67))
%    xh(:,i+1) = A8t*xh(:,i) + B8t*p_vec(:,i);
%    yh(:,i) = C8t*xh(:,i) + D8t*p_vec(:,i);  
%     if(i>app_end) iter8=iter8+1; end
%   elseif ((yh_m>67) && (yh_m<72))
%    xh(:,i+1) = A9t*xh(:,i) + B9t*p_vec(:,i);
%    yh(:,i) = C9t*xh(:,i) + D9t*p_vec(:,i);  
%     if(i>app_end) iter9=iter9+1; end 
%   elseif ((yh_m>72) &&(yh_m<77))
%    xh(:,i+1) = A10t*xh(:,i) + B10t*p_vec(:,i);
%    yh(:,i) = C10t*xh(:,i) + D10t*p_vec(:,i);  
%     if(i>app_end) iter10=iter10+1; end    
%   elseif (yh_m>77)
%    xh(:,i+1) = A11t*xh(:,i) + B11t*p_vec(:,i);
%    yh(:,i) = C11t*xh(:,i) + D11t*p_vec(:,i);  
%     if(i>app_end) iter11=iter11+1; end  
%   end
  
  
 end
toc
%app_end=1;
%  err_diff = abs(y_temp_vec(:,app_end:ver_end)-yh(:,app_end:ver_end) - offset)./(y_temp_vec(app_end:ver_end)-273.15);
err_diff = abs(y_temp_vec(app_end:ver_end,:)'-yh(:,app_end:ver_end) - offset)./(y_temp_vec(app_end:ver_end,:)'-273.15);%-273.15);
  mean_err = mean(err_diff');
  max_err = max(mean_err);
% yh(1) = C1*xh(:,1) + D1*p_vec(1);


plot(t(1:ver_end)',yh(1,1:ver_end)+20,'--r','linewidth',2); hold on;


legend('reference','sid');
xlabel('time (sec)');
ylabel('temp (celcius)');
title('MIMO response');


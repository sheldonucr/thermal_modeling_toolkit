function [data_er,data_er2] = cross_valid_plot(yy0,y_ssm1,y_ssmm1,time)

figure;
plot(time,yy0(:,1),'k','linewidth',2); hold on;
plot(time,y_ssm1{1,1}.OutputData(:,1),':r','linewidth',2); hold on;
plot(time,y_ssmm1{1,1}.OutputData(:,1),'--g','linewidth',2); hold on;
%plot(times,y_ssmm1(:,1),'--g','linewidth',2); hold on;
%legend('verif','sid','thermsid');
legend('train','sid','thermsid');
xlabel('time (sec)');
ylabel('temp (celcius)');
title('Temperature Response Core 0');

figure;
plot(time,yy0(:,2),'k','linewidth',2); hold on;
plot(time,y_ssm1{1,1}.OutputData(:,2),':r','linewidth',2); hold on;
plot(time,y_ssmm1{1,1}.OutputData(:,2),'--g','linewidth',2); hold on;
%plot(times,y_ssmm1(:,2),'--g','linewidth',2); hold on;
%legend('verif','sid','thermsid');
legend('train','sid','thermsid');
xlabel('time (sec)');
ylabel('temp (celcius)');
title('Temperature Response Core 1');

figure;
plot(time,yy0(:,3),'k','linewidth',2); hold on;
plot(time,y_ssm1{1,1}.OutputData(:,3),':r','linewidth',2); hold on;
plot(time,y_ssmm1{1,1}.OutputData(:,3),'--g','linewidth',2); hold on;
%plot(times,y_ssmm1(:,3),'--g','linewidth',2); hold on;
%legend('verif','sid','thermsid');
legend('train','sid','thermsid');
xlabel('time (sec)');
ylabel('temp (celcius)');
title('Temperature Response Core 2');

figure;
plot(time,yy0(:,4),'k','linewidth',2); hold on;
plot(time,y_ssm1{1,1}.OutputData(:,4),':r','linewidth',2); hold on;
plot(time,y_ssmm1{1,1}.OutputData(:,4),'--g','linewidth',2); hold on;
%plot(times,y_ssmm1(:,4),'--g','linewidth',2); hold on;
%legend('verif','sid','thermsid');
legend('train','sid','thermsid');
xlabel('time (sec)');
ylabel('temp (celcius)');
title('Temperature Response Core 3');

figure;
plot(time,yy0(:,5),'k','linewidth',2); hold on;
plot(time,y_ssm1{1,1}.OutputData(:,5),':r','linewidth',2); hold on;
plot(time,y_ssmm1{1,1}.OutputData(:,5),'--g','linewidth',2); hold on;
%plot(times,y_ssmm1(:,5),'--g','linewidth',2); hold on;
%legend('verif','sid','thermsid');
legend('train','sid','thermsid');
xlabel('time (sec)');
ylabel('temp (celcius)');
title('Temperature Response Cache');

ya = y_ssm1{1,1}.OutputData(:,:);
yb = y_ssmm1{1,1}.OutputData(:,:);

sd_ssm = abs( yy0 - y_ssm1{1,1}.OutputData(:,:) );
sd_ssmm = abs( yy0 - y_ssmm1{1,1}.OutputData(:,:) );
%sd_ssmm = abs( yy0 - y_ssmm1 );
pc_ssm = [];
pc_ssmm = [];
for j = 1:length(yy0)
    pc_ssm = [ pc_ssm; [ sd_ssm(j,1)/yy0(j,1) sd_ssm(j,2)/yy0(j,2) sd_ssm(j,3)/yy0(j,3) sd_ssm(j,4)/yy0(j,4) sd_ssm(j,5)/yy0(j,5) ] ];
    pc_ssmm = [ pc_ssmm; [ sd_ssmm(j,1)/yy0(j,1) sd_ssmm(j,2)/yy0(j,2) sd_ssmm(j,3)/yy0(j,3) sd_ssmm(j,4)/yy0(j,4) sd_ssmm(j,5)/yy0(j,5) ] ];
end

data_er = [];
data_er2 = [];
for j = 1:size(pc_ssm,2)
    data_er = [ data_er; [ max(pc_ssm(:,j)) max(pc_ssmm(:,j)) mean(pc_ssm(:,j)) mean(pc_ssmm(:,j)) std(pc_ssm(:,j)) std(pc_ssmm(:,j)) ] ];
    data_er2 = [ data_er2; [ max(sd_ssm(:,j)) max(sd_ssmm(:,j)) mean(sd_ssm(:,j)) mean(sd_ssmm(:,j)) std(sd_ssm(:,j)) std(sd_ssmm(:,j)) ] ];
end

data_er = 100*data_er;
data_er2 = 100*data_er2;
    
    
% figure;
% plot(time,uu(:,1),'k');hold on;
% plot(time,yy(:,1),'b');hold off;
% 
% figure;
% plot(times,uu0(:,1),'k');hold on;
% plot(times,yy0(:,1),'k');hold off;


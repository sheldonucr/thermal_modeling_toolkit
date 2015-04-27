close all
clear all

H = load('all_data.txt');

for i = 1:length(H)/6
    times(i)= H(6*(i-1)+1) - round(H(1));
    val1(i) = H(6*(i-1)+2); % core0 temp
    val2(i) = H(6*(i-1)+3); % core1 temp
    val3(i) = H(6*(i-1)+4); % core2 temp
    val4(i) = H(6*(i-1)+5); % core3 temp
    val5(i) = H(6*(i-1)+6); % cache temp
end

load('pr_core0.mat');

start1 = start(1);
start2 = start(2);
start3 = start(3);
start4 = start(4);
start5 = start(5);

% pr_data = load('poles_residues.mat')

times1 = times(start(1):length(times));
times2 = times(start(2):length(times));
times3 = times(start(3):length(times));
times4 = times(start(4):length(times));
times5 = times(start(5):length(times));
v1 = val1(start1:length(val1));
v2 = val1(start2:length(val1));
v3 = val1(start3:length(val1));
v4 = val1(start4:length(val1));
v5 = val1(start5:length(val1));

[step_res1, orig_step_res1] = compute_step_res(times1, v1, si1{1}, b1{1}, L(1));
[step_res2, orig_step_res2] = compute_step_res(times2, v2, si1{2}, b1{2}, L(2));
[step_res3, orig_step_res3] = compute_step_res(times3, v3, si1{3}, b1{3}, L(3));
[step_res4, orig_step_res4] = compute_step_res(times4, v4, si1{4}, b1{4}, L(4));
[step_res5, orig_step_res5] = compute_step_res(times5, v5, si1{5}, b1{5}, L(5));

orig_step_res1 = [zeros(1,start1-1), orig_step_res1];
orig_step_res2 = [zeros(1,start2-1), orig_step_res2];
orig_step_res3 = [zeros(1,start3-1), orig_step_res3];
orig_step_res4 = [zeros(1,start4-1), orig_step_res4];
orig_step_res5 = [zeros(1,start5-1), orig_step_res5];

% N = 100
% no = length(times);
% time_intval = times(no) - times(1);
% time_step = time_intval / no;
% % t = 0:time_step:times(no)
% for i = 1:no
%     t(i) = (i-1)*time_step
% end
% 
% orig_step_res1 = zeros(1, no);
% orig_step_res2 = zeros(1, no);
% orig_step_res3 = zeros(1, no);
% orig_step_res4 = zeros(1, no);
% orig_step_res5 = zeros(1, no);
% for i = 1 : length(core0.si1)
%     orig_step_res1 = orig_step_res1 + 20*(core0.b1(i)/core0.si1(i)) * (exp(core0.si1(i) * t) - 1);
%     orig_step_res2 = orig_step_res2 + 20*(core1.b1(i)/core1.si1(i)) * (exp(core1.si1(i) * t) - 1);
%     orig_step_res3 = orig_step_res3 + 20*(core2.b1(i)/core2.si1(i)) * (exp(core2.si1(i) * t) - 1);
%     orig_step_res4 = orig_step_res4 + 20*(core3.b1(i)/core3.si1(i)) * (exp(core3.si1(i) * t) - 1);
%     orig_step_res5 = orig_step_res5 + 20*(cache.b1(i)/cache.si1(i)) * (exp(cache.si1(i) * t) - 1);
% end;

% step_res_all = step_res1 + step_res2 + step_res3 + step_res4 + step_res5;
orig_step_res_all = val1(1) + orig_step_res1 + orig_step_res2 + orig_step_res3 + orig_step_res4 + orig_step_res5;

% plot_waveform(times, val1, step_res_all, orig_step_res_all, L(1));

figure
hold on
plot(times, val1,'r');
plot(times, real(orig_step_res_all), 'b:');
legend('original','computed model');
hold off


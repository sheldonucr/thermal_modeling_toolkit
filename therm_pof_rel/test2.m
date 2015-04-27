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

core0 = load('poles_residues0.mat')
core1 = load('poles_residues1.mat')
core2 = load('poles_residues2.mat')
core3 = load('poles_residues3.mat')
core4 = load('poles_residues4.mat')

% L = 100
% [step_res1, orig_step_res1] = compute_step_res(times, val1, core0.si1, core0.b1, L)
% [step_res2, orig_step_res2] = compute_step_res(times, val2, core1.si1, core1.b1, L)
% [step_res3, orig_step_res3] = compute_step_res(times, val3, core2.si1, core2.b1, L)
% [step_res4, orig_step_res4] = compute_step_res(times, val4, core3.si1, core3.b1, L)
% [step_res5, orig_step_res5] = compute_step_res(times, val5, core4.si1, core4.b1, L)

no = length(times);
t_offset = -log(times(1));
for i = 1 : no
    t(i) = log(times(i)) + t_offset;
end
orig_step_res1 = zeros(1, no) + 35;
orig_step_res2 = zeros(1, no) + 35;
orig_step_res3 = zeros(1, no) + 35;
orig_step_res4 = zeros(1, no) + 35;
orig_step_res5 = zeros(1, no) + 35;
for i = 1 : length(core0.si1)
    orig_step_res1 = orig_step_res1 + 20*(core0.b1(i)/core0.si1(i)) * (exp(core0.si1(i) * t) - 1);
end
for i = 1 : length(core1.si1)
    orig_step_res2 = orig_step_res2 + 20*(core1.b1(i)/core1.si1(i)) * (exp(core1.si1(i) * t) - 1);
end
for i = 1 : length(core2.si1)
    orig_step_res3 = orig_step_res3 + 20*(core2.b1(i)/core2.si1(i)) * (exp(core2.si1(i) * t) - 1);
end
for i = 1 : length(core3.si1)
    orig_step_res4 = orig_step_res4 + 20*(core3.b1(i)/core3.si1(i)) * (exp(core3.si1(i) * t) - 1);
end
for i = 1 : length(core4.si1)
    orig_step_res5 = orig_step_res5 + 20*(core4.b1(i)/core4.si1(i)) * (exp(core4.si1(i) * t) - 1);
end;

orig_step_res0 = orig_step_res1 + orig_step_res2 + orig_step_res3 + orig_step_res4 + orig_step_res5

figure
hold on
plot(times, val1,'r');
plot(times, orig_step_res0 - 35*4, 'b:');
legend('original','computed model');
hold off


close all
clear all

for i = 1:5
    inputfile(i,:) = ['core', int2str(i-1), '_data.txt'];
    H(i,:) = load(inputfile(i,:), '-ascii');
end

times = cell(3);
val1 = cell(3);
val2 = cell(3);
val3 = cell(3);
val4 = cell(3);
val5 = cell(3);
for j = 1:5
    for i = 1:length(H(j, :))/6
        times{j}(i) = H(j, 6*(i-1)+1) - round(H(j, 1)); % time
        val1{j}(i) = H(j, 6*(i-1)+2); % core0 temp
        val2{j}(i) = H(j, 6*(i-1)+3); % core1 temp
        val3{j}(i) = H(j, 6*(i-1)+4); % core2 temp
        val4{j}(i) = H(j, 6*(i-1)+5); % core3 temp
        val5{j}(i) = H(j, 6*(i-1)+6); % cache temp
    end
    len(j) = length(H(j, :))/6;
end

M = 50;
% L = [100 100 100 100 100];
L(1) = 100;
L(2) = 78;
L(3) = 84;
L(4) = 88;
L(5) = 70;

start(1) = 1;
start(2) = 1;
start(3) = 1;
start(4) = 1;
start(5) = 1;

% start(1) = 1;
% start(2) = 4;
% start(3) = 5;
% start(4) = 5;
% start(5) = 6;

t= cell(3);
v1 = cell(3);
v2 = cell(3);
v3 = cell(3);
v4 = cell(3);
v5 = cell(3);

for i = 1:5
    t{i} = times{i}(start(i):len(i));
    v1{i} = val1{i}(start(i):len(i));
    v2{i} = val2{i}(start(i):len(i));
    v3{i} = val3{i}(start(i):len(i));
    v4{i} = val4{i}(start(i):len(i));
    v5{i} = val5{i}(start(i):len(i));
end

si1 = cell(3);
b1 = cell(3);
[si1{1}, b1{1}] = matpencil_intel_run(t{1}, v1{1}, M, L(1));
[si1{2}, b1{2}] = matpencil_intel_run(t{2}, v1{2}, M, L(2));
[si1{3}, b1{3}] = matpencil_intel_run(t{3}, v1{3}, M, L(3));
[si1{4}, b1{4}] = matpencil_intel_run(t{4}, v1{4}, M, L(4));
[si1{5}, b1{5}] = matpencil_intel_run(t{5}, v1{5}, M, L(5));

for i = 1:5
    offset(i) = times{i}(start(i));
end

outputfile = ['pr_core0','.mat'];
save (outputfile, 'M', 'L', 'si*', 'b*','offset','start');

% for checking the poles and the residues

step_res = cell(3);
orig_step_res = cell(3);
for i = 1:5
    [step_res1{i}, orig_step_res1{i}] = compute_step_res(t{i}, v1{i}, si1{i}, b1{i}, L(i));
end

% extend val1
% for i = 1:5
%     val1{i} = [zeros(1,start{i}-1) + 35, val1{i}];
% end

% only for debugging
% step_res1 = [zeros(1,start1-1), step_res1, zeros(1,start5-start1) + step_res1(length(step_res1))]
% step_res2 = [zeros(1,start2-1), step_res2, zeros(1,start5-start2) + step_res2(length(step_res2))]
% step_res3 = [zeros(1,start3-1), step_res3, zeros(1,start5-start3) + step_res3(length(step_res3))]
% step_res4 = [zeros(1,start4-1), step_res4, zeros(1,start5-start4) + step_res4(length(step_res4))]
% step_res5 = [zeros(1,start5-1), step_res5, zeros(1,start5-start5) + step_res5(length(step_res5))]

% for i = 1:5
%     step_res1{i} = [zeros(1,start(i)-1), step_res1{i}, zeros(1,start(5)-start(i)) + step_res1{i}(length(step_res1{i}))];
% end

% orig_step_res1 = [zeros(1,start1-1) + 35, orig_step_res1]
% orig_step_res2 = [zeros(1,start2-1) + 35, orig_step_res2]
% orig_step_res3 = [zeros(1,start3-1) + 35, orig_step_res3]
% orig_step_res4 = [zeros(1,start4-1) + 35, orig_step_res4]
% orig_step_res5 = [zeros(1,start5-1) + 35, orig_step_res5]
% for i = 1:5
%     orig_step_res1{i} = [zeros(1,start(i)-1) + 35, orig_step_res1{i}];
% end

for i = 1:5
    plot_waveform(t{i}, v1{i}, step_res1{i}, orig_step_res1{i} + v1{i}(1), L(i));
end 
% figure
% hold on
% plot(times(1:length(times)-1), val0(1:length(val0)-1),'r');
% plot(times(1:length(times)-1), orig_step_res0(1:length(orig_step_res0)-1), 'b:');
% legend('original','computed model');
% hold off

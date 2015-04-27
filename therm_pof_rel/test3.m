clear all
close all

H = load('random_data.txt');
for i = 1:size(H,1)
    times(i)= H(i,1) - round(H(1,1));
    val1(i) = H(i,2); % core0 temp
    val2(i) = H(i,3); % core1 temp
    val3(i) = H(i,4); % core2 temp
    val4(i) = H(i,5); % core3 temp
    val5(i) = H(i,6); % cache temp
end

nSegments = 10;
nCores = 5;
tLength = [100,100,101,101,101,101,101,101,101,101];
tIndex = zeros(1,11);
tIndex(1) = 1;
for i = 1:nSegments
    tIndex(i+1) = tIndex(i) + tLength(i);
end

P = zeros(nSegments+1,5);
for i = 1:nSegments
    for j = 1:nCores
        P(i,j) = H(tIndex(i),j+6)/20;
        if (P(i,j) == 0)
            P(i,j) = -1;
        end
    end
end

plot(times, val1);

load('pr_core0.mat');

nCell = ceil(nSegments^(1/2));
t = cell(nCell);
v1 = cell(nCell);
v2 = cell(nCell);
v3 = cell(nCell);
v4 = cell(nCell);
v5 = cell(nCell);

for i = 1:nSegments
    t{i} = times(tIndex(i):tIndex(i+1)-1) - (i-1)*0.1; 
    t{i}(1) = round(t{i}(1));
    if (i == 2 | i == 3)
        t{i} = t{i} + 0.001;
        t{i} = [0, t{i}];
    elseif (i > 3)
        t{i} = t{i} + 0.000991;
        t{i} = [0, t{i}];
    end
    v1{i} = val1(tIndex(i):tIndex(i+1)-1);
    v2{i} = val2(tIndex(i):tIndex(i+1)-1);
    v3{i} = val3(tIndex(i):tIndex(i+1)-1);
    v4{i} = val4(tIndex(i):tIndex(i+1)-1);
    v5{i} = val5(tIndex(i):tIndex(i+1)-1);
end

for i = 1:nCores
    mindist = 100;
    for j = 1:10
        if (abs(offset(i) - times(j)) < mindist)
            start(i) = j;
            mindist = abs(offset(i) - times(j));
        end
    end
end

step_res1 = cell(nSegments,nCell);
orig_step_res1 = cell(nSegments,nCell);

for i = 1:nSegments
    for j = 1:nCores
        ti = t{i}(start(j):length(t{i}));
        ti(1) = offset(j);
        v1i = v1{i}(start(j):tLength(i));
        [step_res1{i,j}, orig_step_res1{i,j}] = compute_step_res4(ti, v1i, si1{j}, b1{j}, L(j), P(:,j), 20, i);
    end
end

step_res1_all = cell(nCell);
orig_step_res1_all = cell(nCell);

for i = 1:nSegments
    for j = 1:nCores
        orig_step_res1{i,j} = [zeros(1,start(j)-1)+orig_step_res1{i,j}(1), orig_step_res1{i,j}];
    end
    orig_step_res1_all{i} = 35 + orig_step_res1{i,1} + orig_step_res1{i,2} + orig_step_res1{i,3} + orig_step_res1{i,4} + orig_step_res1{i,5};
end

for i = 1:10
%     figure
%     hold on
    if (i > 1)
        t{i} = t{i}(2:length(t{i}));
        orig_step_res1_all{i} = orig_step_res1_all{i}(2:length(orig_step_res1_all{i}));
    end   
%     plot(t{i}, v1{i},'r');
%     plot(t{i}, real(orig_step_res1_all{i}), 'b:');
%     legend('original','computed model');
%     hold off
end

v1_all = [];
v1_res = [];
for i = 1:10
    v1_all = [v1_all, v1{i}];
    v1_res = [v1_res, orig_step_res1_all{i}];
end
times = times(1:length(times)-1);
figure
hold on
plot(times, v1_all,'r');
plot(times, real(v1_res), 'b:');
legend('original','computed model');
hold off

% v1_diff = v1_all - real(v1_res);
% figure
% hold on
% plot(times, v1_diff,'r');
% legend('difference');
% hold off

% % P1 = P;
% P(11,:) = 1;
% for j = 1:5
%     P1(1,j) = P(1,j);
%     P1(2,j) = P(1,j);
%     for i = 2:11
%         P1(2*i-1,j) = P(i-1,j);
%         P1(2*i,j) = P(i,j);
%     end
% end
% P1 = 10*P1;
% for i = 1:5
%     P1(:,i) = P1(:,i) + i*40;
% end
% 
% pt(1) = 0;
% pt(2) = 0;
% for i = 2:11
%     pt(2*i-1) = 0.1*(i-1) - 1e-100;
%     pt(2*i) = 0.1*(i-1);
% end
% figure
% hold on 
% for i = 1:5
%     plot(pt,P1(1:22,i), 'b');
% end
% % legend('Power input');
% hold off


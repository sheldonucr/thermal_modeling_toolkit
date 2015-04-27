function plot_waveform(times, val, step_res, orig_step_res, L)

% M = 50;
% L = 80;
N = 2*L+1;
no = length(val);
t_offset = -log(times(1));
for i=1:no
    times_log(i) = log(times(i)) + t_offset;
end
time_intval = times_log(no) - times_log(1);
time_step = time_intval / (N-1);
t = 0:time_step:times_log(no);

figure;
plot(times_log, val-val(1), 'r'); 
hold on; 
plot(t, real(step_res), 'b:');
legend('original','computed');

% compute the response in the original time scale
figure
hold on
plot(times(1:length(times)-1), val(1:length(val)-1),'r');
plot(times(1:length(times)-1), real(orig_step_res(1:length(orig_step_res)-1)), 'b:');
legend('original','computed model');
hold off


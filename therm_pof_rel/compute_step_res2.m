function [step_res, orig_step_res] = compute_step_res2(times, val, si, b, L)

% M = 50;
% L = 80;
N = 2*L+1;
no = length(val);

org_times = times;
org_val = val;

t_offset = -log(times(1));
for i=1:no
    times(i) = log(times(i)) + t_offset;
end

time_intval = times(no) - times(1);
time_step = time_intval / (N-1);
t = 0:time_step:times(no);

% finally we compare the step responses from the new transfer function.
step_res = zeros(1, length(t));
for i = 1 : length(si)
    step_res = step_res + 20*(b(i)/si(i)) * (exp(si(i) * t) - 1);
end;

% compute the response in the original time scale
for i = 1 : no
    scaled_t(i) = log(org_times(i)) + t_offset;
end
orig_step_res = zeros(1, no);
for i = 1 : length(si)
    orig_step_res = orig_step_res + 20*(b(i)/si(i)) * (exp(si(i) * scaled_t) - 1);
end;
% orig_step_res = val(1) + orig_step_res;
orig_step_res = orig_step_res;

% function [si, b, impulse_series, step_res, orig_step_res] = matpencil_intel_run(times, val)
function [si, b] = matpencil_intel_run(times, val, M, L)

% In the function, we assume the given data are step responses, as a
% result, we need to compute the impluse responses from the step response
% as matrix pencil assumes that the given data is the impulse responses.

%
% * $RCSfile: matpencil_intel_run.m,v $
% * $Revision: 1.1.1.1 $
% * $Date: 2006/12/08 23:14:28 $, (c) Copyright 2006- by University of California
% * Authors: Sheldon Tan
% 

% This algorithm is implemented based on the following paper:
% Y. Hua and T. Sarkar, "Generalized pencil-of-function method for
% extracting poles of an EM system from its transient response", IEEE
% Trans. on Antennas and Propagation, vol. 37, no.2 pp. 229-234, 1989.

% Some parameters for the matrix pencil call

% L: window size
% M: number of poles used to model the input
% N: total number of samplings

% hold on;

% M = 50;
% L = 100;
N = 2*L+1;

% one example

% H = load(file_name);
% H = load('core1_data.txt');
% %H = H(1:600,1);
% 
% % The data comes from the winspice "write file variable" command output.
% % The output needs to be preprocessed so that the indices for each time
% % step are removed.
% 

%%%%%%%%%%%%% Sampling in time domain %%%%%%%%%%%%%%%%
no = length(val);

org_times = times;
org_val = val;

t_offset = -log(times(1));
for i=1:no
    times(i) = log(times(i)) + t_offset;
end


time_intval = times(no) - times(1);
time_step = time_intval / N;
t = 0:time_step:times(no);

disp(sprintf('number of time steps in the original data: %d\n', no));
disp(sprintf('number of time steps used for sampling: %d\n', N));


% create different time scale to verify the fitting results.
time_step1 = 0.5*time_intval/N;
t1 = 0:time_step1:2*times(no);

for i=1:no
    val(i) = val(i) - 35;
end

% plot(times, val,'b');

% fit the required sampling points based on the given pwl input
for i=1:N
    time_series(i) = times(1) + ((i-1)*time_step);
    tm = 1;
    while(tm <= no)
        if(time_series(i) >= times(tm))
            start = tm;
        end;
        tm = tm + 1;
    end;
    
    tm = 1;
    while(tm <= no)
        if(time_series(i) <= times(no-tm+1))
            endd = no - tm + 1;
        end;
        tm = tm + 1;
    end;
   
    if(start == endd)
        amp_series(i) = val(start);
    else
        amp_series(i) = val(start) + (((val(endd)-val(start))/(times(endd)-times(start)))*(time_series(i)-times(start)));
    end;
    
end;

% plot(time_series,amp_series,'g:+');

%%%%%%%%%%%%%%%%%%%%%%%%% compute the impluse responses from the step responses

% we just take the numerical derivative of the given data.
impulse_series(1) = 0;
for i=2:N
    impulse_series(i) = 1/20*(amp_series(i) - amp_series(i-1))/time_step;
end

% figure
% plot(time_series, impulse_series,'b');
% hold on
% pp = spline(time_series, hs);
% 
% new_N = 1000;
% time_intval = times(no) - times(1);
% new_time_step = time_intval / new_N;
% new_time = 0:new_time_step:times(no);
% 
% % the smoothed values
% smooth_serise = ppval(pp, new_time);

% disp(sprintf('number of time steps used for smoothing: %d\n', new_N));

% now, we compute the integration area
area_1 = 0;
for i=1:N-1
    area_1 = area_1 + (impulse_series(i)+impulse_series(i+1))*time_step/2;
end
% area_2 = 0;
% for i=1:new_N-1
%     area_2 = area_2 + (smooth_series(i)+smooth_series(i+1))*new_time_step/2;
% end

disp(sprintf('The integration area for un-smooth impluse curve is %f\n',area_1));
% disp(sprintf('this integration area for smoothed impluse curve is %f\n',area_2));


%%%%%%%%%%% Y matrix construction %%%%%%%%%%%%%%%

Y1 = zeros((N-L),L);
Y2 = zeros((N-L),L);

% build the Y1 matrix
for j=1:L
    for i=1:(N-L)
        %Y1(i,j) = amp_series(j+i-1);
        Y1(i,j) = impulse_series(j+i-1);
    end;
end;

% build the Y2 matrix
for j=1:L
    for i=1:(N-L)
        %Y2(i,j) = amp_series(j+i);
        Y2(i,j) = impulse_series(j+i);
    end;
end;

[U,D,V] = svds(Y1, M);

%[U,D,V] = svd(Y1);

%%%%%%%%%%%% Eigen (Pole) Extraction %%%%%%%%%%%%%
for i = 1:M
    D_inv(i, i) = diag(1/D(i,i));
end;

Z =  D_inv * conj(U)' * Y2 * V;

z_poles = conj(eig(Z))';

%%%%%%%%%%%%%%% Calculate for residue %%%%%%%%%%%%%%%%%%%%%
Z0 = diag(z_poles);

% build the Z1 matrix
for i = 0: N-L-1
    for j = 1: M
        Z1(i + 1, j) = z_poles(j) ^ i;
    end;
end;

% build the Z2 matrix
for i = 0: L - 1
    for j = 1: M
        Z2(j, i+1) = z_poles(j) ^ i;
    end;    
    %C2(:, i+1) = (f.').^i;
end;

%  Y1 = Z1*B*Z2 (eq.5)
%  Y2 = Z1*B*Z0*Z2 (eq.6)
%  so b should equal to b2 in the following
B = (Z1\Y1)/Z2;
% B2 = (Z1\Y2)/(Z0*Z2);
b = diag(B);
% b2 = diag(B2);

% compute the time domain poles (si) which is related to Z-domain poles by
% zi = exp(si*time_step);
for i = 1:M
    si(i) = log(z_poles(i))/time_step;
end;

tmp = zeros(1, length(t));
for i = 1 : length(si)
    tmp = tmp + b(i) * exp(si(i) * t);
end;
% plot(t, tmp, 'g:');

tmp1 = zeros(1, length(t1));
for i = 1 : length(si)
    tmp1 = tmp1 + b(i) * exp(si(i) * t1);
end;
% plot(t1, tmp1, 'b:');
% legend('original','fitting time scale','another time scale');
% hold off

% finally we compare the step responses from the new transfer function.
% figure;
% plot(times, val,'r'); 
% hold on; 

% for i=1:length(si)
%     step_res(i) = step_res(i) + 35;
% end
% 

step_res = zeros(1, length(t));
for i = 1 : length(si)
    step_res = step_res + 20*(b(i)/si(i)) * (exp(si(i) * t) - 1);
end;
% plot(t, step_res, 'g:');
% legend('original','computed');

% compute the response in the original time scale
for i = 1 : no
    scaled_t(i) = log(org_times(i)) + t_offset;
end

orig_step_res = zeros(1, no);
for i = 1 : length(si)
    orig_step_res = orig_step_res + 20*(b(i)/si(i)) * (exp(si(i) * scaled_t) - 1);
end;

orig_step_res = 35 + orig_step_res;

% figure
% hold on
% plot(org_times(1:no-1), org_val(1:no-1),'b');
% plot(org_times(1:no-1), orig_step_res(1:no-1), 'g:+');
% 
% legend('original','computed model');
% hold off

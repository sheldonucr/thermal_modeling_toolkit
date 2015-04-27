function  [G, C, B, L, U] = rlz_vf3(AA, BB, CC, DD, inNum, outNum, file_name, subckt_name, is_org, is_sub_ckt, is_new_file, sim_opts)


N=size(AA,1);
cindex=zeros(1,N);
for m=1:N-1 
  if (AA(m,m+1))~=0  
    if m==1 
      cindex(m)=1;
    else
      if cindex(m-1)==0 || cindex(m-1)==2
        cindex(m)=1; cindex(m+1)=2; 
      else
        cindex(m)=2;
      end
    end 
  end
end
p2=[];
for m=1:N
    if(cindex(m) ~=0)
        p2(m)=(-1)^(cindex(m)-1);
    else
        p2(m)=1;
    end
end
p2 = sparse(diag(p2));

G = p2 * AA;
C = -p2 * eye(size(AA));
B = -p2 * BB;
L = CC;
D = DD;

%  G = AA;
%  C = - eye(size(AA));
%  B = - BB;
%  L = CC;
%  D = DD;


%scale = 1e10;
scale =-1e6 ;

G = G./scale;
C = C./scale;
B = B./scale;
% for i = 1:6
%     I_face{i} = I_face{i}./scale;
% end
% some information about the reduced circuit.
n_dim = size(B,1);

if (is_org == 0)
    % Diagonalize G, C
    [U, D] = eig(C, G);
    C = U'*C*U;
    G = U'*G*U;
    B = U'*B;
    L = L*U;
    % zero the small number caused by numerical error
    for i = 1:n_dim
        for j = 1:n_dim
            if (i ~= j)
                C(i,j) = 0;
                G(i,j) = 0;
            end
        end
    end
    % U = [];
else
    U = [];
end

% output
if is_sub_ckt == 0
    [fid, msg] = fopen(file_name,'w');
else
    if is_new_file == 1
        [fid, msg] = fopen(file_name,'w');
    else
        [fid, msg] = fopen(file_name,'a');
    end
end

if (fid == -1)
    printf('Fatal error: $s\n',msg);
    return;
end

if is_org == 1
    fprintf(fid,'** Original thermal ckt realization\n');
    if is_sub_ckt == 1
        fprintf(fid,'.SUBCKT %s ', subckt_name);
    end
else
    fprintf(fid,'** Reduced thermal ckt realization\n');
    if is_sub_ckt == 1
        fprintf(fid,'.SUBCKT %s ', subckt_name);
    end
end
port_name = cell(inNum,1);
for i = 1:inNum
    port_name{i} = strcat('port_',num2str(i));
    if is_sub_ckt == 1
        if mod(i,10) == 0 && i ~= inNum
            fprintf(fid, '%s\n+ ', port_name{i});
        else
            fprintf(fid, '%s', port_name{i}, ' ');
        end
    end
end
fprintf(fid,'\n');

% capacitor realization
k = 1;
for i = 1:n_dim
    for j = i+1:n_dim
        if  C(i,j) ~= 0
            Cap_s(k) = i;
            Cap_e(k) = j;
            Cap(k) = -C(i,j);
            Cap_name{k} = strcat('C',num2str(k));
            k = k+1;
        end        
    end       
end

for i = 1:n_dim
    sum = 0;
    for j = 1:n_dim
        sum = sum + C(i,j);
    end
    if sum ~= 0
        Cap_s(k) = i;
        Cap_e(k) = 0;
        Cap(k) = sum;
        Cap_name{k} = strcat('C',num2str(k));
        k = k+1;
    end
end

NC = k - 1;
for k=1:NC
    fprintf(fid,'%s %d %d %2.15e\n',Cap_name{k},Cap_s(k),Cap_e(k),Cap(k));
end

% resistor realization
k = 1;

for i = 1:n_dim
    for j = i+1:n_dim
        if  G(i,j) ~= 0
            R_s(k) = i;
            R_e(k) = j;
            R(k) = -1/G(i,j);    
            R_name{k} = strcat('R',num2str(k));
            k = k+1;
        end        
    end       
end

for i = 1:n_dim
    sum = 0;
    for j = 1:n_dim
        sum = sum + G(i,j);
    end
    if sum ~= 0
        R_s(k) = i;
        R_e(k) = 0;
        R(k) = 1/sum;
        R_name{k} = strcat('R',num2str(k));
        k = k+1;
    end
end

NR = k - 1;
for k=1:NR
    fprintf(fid,'%s %d %d %2.15e\n',R_name{k},R_s(k),R_e(k),R(k));
end

% some additional Resistance, for output voltage 
for i = 1:inNum
    V_ex_e{i} = port_name{i};
   % V_ex_s{i} = strcat('n_0_',num2str(i)); 
   %V_ex_e is connect to input current source, while V_ex_s is connect to
   %output voltage source
    V_ex_s{i} = strcat('o_port_', num2str(i));
    V_ex_name{i} = strcat('V_ex_',num2str(i));
end
% 
if outNum >inNum
  for i=inNum+1:outNum
    %V_ex_e{i} = 0;
    V_ex_s{i} = strcat('o_port_', num2str(i));
   % V_ex_name{i} = strcat('V_ex_',num2str(i));
  end
end
%%%if inNum>OutNum, just ignore that output


%for i=1:max(outNum,inNum)
for i=1:inNum
    fprintf(fid,'%s %s %s 0\n',char(V_ex_name{i}),char(V_ex_s{i}),char(V_ex_e{i}));
end 

% CCCS realization
k = 1;
for i = 1:n_dim
    for j = 1:inNum
        if (B(i,j) ~= 0)
            F_s(k) = i;
            F_e(k) = 0;
            F_elem{k} = V_ex_name{j};
            F_val(k) = B(i,j);
            F_name{k} = strcat('F_',num2str(i),'_',num2str(j));
            k = k+1;
        end
    end
end
n_cccs = k-1;
for i=1:n_cccs
    fprintf(fid,'%s %d %d %s %2.15e\n',char(F_name{i}),F_s(i),F_e(i),F_elem{i},F_val(i));        
end

% VCVS realization
k = 1;
if nnz(D) == 0
    is_zero_D = 1;
else
    is_zero_D = 0;
end
kk = 0;
for j = 1:outNum
   % col_nnz = nnz(B(:,j));
    row_nnz = nnz(L(j,:));
    nz_idx = 0;
    for i = 1:n_dim        
        if (L(j,i) ~= 0)%(B(i,j) ~= 0)
            nz_idx = nz_idx+1;
            if nz_idx == 1
                E_s{k} = V_ex_s{j}; %connect to input via V_ex
            else
                E_s{k} = E_e{k-1};
            end
            if ( (nz_idx == row_nnz) && (is_zero_D == 1))
                E_e{k} = '0';
            elseif ( (nz_idx == row_nnz) && (is_zero_D == 0))
                    E_e{k} = strcat('n_', num2str(i), '_', num2str(j));
                    kk = kk +1;
                    E_end{kk} = E_e{k};
            else
                   E_e{k} = strcat('n_', num2str(i), '_', num2str(j));
           end
            E_np(k) = i;
            E_nn(k) = 0;
            E_val(k) = L(j,i);
            E_name{k} = strcat('E_',num2str(i),'_',num2str(j));
            k = k+1;
        end
    end
end
n_vcvs = k-1;
for i=1:n_vcvs
    fprintf(fid,'%s %s %s %d %d %2.15e\n',char(E_name{i}),char(E_s{i}),char(E_e{i}),E_np(i),E_nn(i),E_val(i));
end


%CCVS realization
k = 1;
kk = 0;
if(is_zero_D == 0)
 for j = 1:outNum
   % col_nnz = nnz(B(:,j));
    row_nnz = nnz(D(j,:));
    nz_idx = 0;
    for i = 1:inNum        
        if (D(j,i) ~= 0)%(B(i,j) ~= 0)
            nz_idx = nz_idx+1;
            if nz_idx == 1
                kk = kk +1;
                H_s{k} = E_end{kk}; 
            else
                H_s{k} = H_e{k-1};
            end
            if nz_idx == row_nnz
                H_e{k} = '0';
            else
                   H_e{k} = strcat('net_', num2str(i), '_', num2str(j));
            end
            H_elem{k}=V_ex_name{i};
            H_val(k) = -D(j,i);
            H_name{k} = strcat('H_',num2str(i),'_',num2str(j));
            k = k+1;
        end
    end
 end
 
 n_vcvs = k-1;
 for i=1:n_vcvs
    fprintf(fid,'%s %s %s %s %2.15e\n',char(H_name{i}),char(H_s{i}),char(H_e{i}),H_elem{i},H_val(i));
 end
end 





if is_sub_ckt == 1
    fprintf(fid,'.ENDS\n');
else
    % current source realization
    for i = 1:inNum
         I_ex_s{i} = '0';
         I_ex_e{i} = port_name{i};
%         if i < n_port-np+1            
%             I_val{i} = num2str(I_face{1}(1));
%         else
%             I_val{i} = 'PULSE(0V 3.2e+1V 0mS 1mS)';
%         end
       if (i==1)
         I_val{i} = 'AC 1 sin(0 1 1K)';
       else
         I_val{i} = 'AC 0 sin(0 1 1K)';
       end
         %I_val{i} = '1';
         I_ex_name{i} = strcat('I_ex',num2str(i));
           
    end
    for i=1:inNum
        fprintf(fid,'%s %s %s %s\n',char(I_ex_name{i}),char(I_ex_s{i}),char(I_ex_e{i}),I_val{i});
    end
    
    fprintf(fid,'.tran 0.1ms 10ms\n');
   % fprintf(fid, '.AC DEC 100 1 100K\n');
     fprintf(fid, '.AC %s  %s  %s %s\n', char(sim_opts.scanType), char(sim_opts.scanNum), char(sim_opts.freqStart), char(sim_opts.freqEnd));
    fprintf(fid,'.options accurate=0 post=1\n');
    fprintf(fid,'.END\n');
end

fclose(fid);
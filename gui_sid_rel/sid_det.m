function [system] = sid_det(y, u, order, method)
% This is the deterministic system identification function. When the order is
% high enough, the identified system should be exactly the same as the
% original one.
% input:
% y: no x ns, no is the output port number, ns is the sample number
% u: ni x ns, ni is the input port number, ns is the sample number
% order: the order of the identified system
% method: choices are 'alg1', 'alg2' and 'stable'. 'alg1' is the algorithm 1
% in the book; 'alg2' is the algorithm 2 in the book; 'stable' is the
% stabilize algorithm in the book
% output:
% system: the identified system written in matlab state space model
    
    [no, ns] = size(y);
    ni = size(u, 1);
    if size(u, 2) ~= ns
        error('the sample numbers for the input and output are different');
    end
    
    % ii = 2*ceil(order/no); %% i in the book, number of block rows in Hankel
    % matrices
    ii = order+1;
    % ii = ceil(order/no); %% i in the book, number of block rows in Hankel matrices
    ij = ns-2*ii+1; %% j in the book
    
    % build all the hankel matrices from input and ouput data
    yp = build_hank(y(:,1:ii+ij-1), ii); % output previous
    yf = build_hank(y(:,ii+1:2*ii+ij-1), ii); % output future
    ya = [yp;yf]; % output all (previous and future)
    ypp = build_hank(y(:,1:ii+ij), ii+1); % output previous +
    yfm = build_hank(y(:,ii+2:2*ii+ij-1), ii-1); % output future -
    up = build_hank(u(:,1:ii+ij-1), ii); % input previous
    uf = build_hank(u(:,ii+1:2*ii+ij-1), ii); % input future
    ua = [up;uf]; % input all (previous and future)
    upp = build_hank(u(:,1:ii+ij), ii+1); % input previous +
    ufm = build_hank(u(:,ii+2:2*ii+ij-1), ii-1); % input future -
    yii = build_hank(y(:,ii+1:ii+ij),1); % Y_(i|i)
    uii = build_hank(u(:,ii+1:ii+ij),1); % U_(i|i)
    wp = [up;yp]; % input and output previous
    wpp = [upp;ypp]; % input and output previous +
    H = [up;uf;yp;yf]; % the whole handel matrix
%   rank([up;uf])    
%   size([up;uf])    
    % rq (also called lq) decomposition (we use qr on the transpose matrix)
    [Q, R] = qr(H',0); R = R'; % H = R*Q'; Q is actually useless, better NOT computed
    
    % decompose R matrix to get R_yf, R_uf
    idx = zeros(1,6);
    idx(1) = 1; idx(2) = idx(1)+ni*ii; idx(3) = idx(2)+ni; idx(4) = idx(3)+ni*(ii-1); idx(5) = ...
             idx(4)+no*ii; idx(6) = idx(5)+no; idx(7) = idx(6)+no*(ii-1);
    R_dec = cell(6,6); % store decomposed (block) R into cell R_dec
    for j = 1:6
        for i = 1:6

            R_dec{i,j} = R(idx(i):idx(i+1)-1,idx(j):idx(j+1)-1);
        end
    end    
    R_yf = cell2mat(R_dec(5:6,1:6));
    R_uf = cell2mat(R_dec(2:3,1:6));
    R_wp = [cell2mat(R_dec(1,1:6));cell2mat(R_dec(4,1:6))];    
    R_yfm = cell2mat(R_dec(6,1:6)); % R_yf-
    R_ufm = cell2mat(R_dec(3,1:6)); % R_uf-
    R_wpp = [cell2mat(R_dec(1:2,1:6));cell2mat(R_dec(4:5,1:6))]; % R_wp+
    
    proj_p = eye(2*(ni+no)*ii)-R_uf'*pinv(R_uf*R_uf')*R_uf; % projection matrix onto uf
    R_yf_p = R_yf*proj_p; % projected R_yf onto uf
    R_wp_p = R_wp*proj_p; % projected R_wp onto uf
    proj_p_p = eye(2*(ni+no)*ii)-R_ufm'*pinv(R_ufm*R_ufm')*R_ufm; % projection matrix onto uf-
    R_yfm_p = R_yfm*proj_p_p; % projected R_yf- onto uf-
    R_wpp_p = R_wpp*proj_p_p; % projected R_wp+ onto uf-
    
    % the oblique projection to obtain O_i and O_{i-1}
    O_i = R_yf_p*pinv(R_wp_p)*wp; % O_i = yf along uf onto wp
    O_im1 = R_yfm_p*pinv(R_wpp_p)*wpp; % O_{i-1} = yfm along ufm onto wpp
    
    % svd 
    WOW = O_i; % here, we assume identity weight matrix W1 and W2
    [U, S] = svd(WOW);
   %  diag(S) 
    % cut to order; order can be specified by user or determined by the
    % singular values S
    U1 = U(:,1:order); 
    S1 = S(1:order, 1:order); 
    
    % determine Gamma and states
    Gamma_i = U1*S1^0.5; % Gamma_i
    Gamma_im1 = Gamma_i(1:end-no,:); % \underline{Gamma_i}, also Gamma_{i-1}
    Gamma_ip1 = Gamma_i(no+1:end,:); % \overline{Gamma_i}, also Gamma_{i+1}
    x_i = pinv(Gamma_i)*O_i; % X_i
    x_ip1 = pinv(Gamma_im1)*O_im1; % X_(i+1)
    %x_i = (Gamma_i)\O_i; % X_i
    %x_ip1 = (Gamma_im1)\O_im1; % X_(i+1)
    
    
    % identify system matrix, for deterministic system, solve the equation
    % directly
    if strcmp(method,'alg1')
        mat_abcd = [x_ip1;yii]*pinv([x_i;uii]);
        %mat_abcd = [x_ip1;yii]/[x_i;uii];
        mat_a = mat_abcd(1:order,1:order); % A matrix
        mat_d = mat_abcd(order+1:end,order+1:end); % D matrix
        mat_b = mat_abcd(1:order,order+1:end); % B matrix
        mat_c = mat_abcd(order+1:end,1:order); % C matrix
    else
        mat_c = Gamma_i(1:no,:);
        [Q_tmp, R_tmp] = qr([Gamma_i,rand(no*ii,no*ii-order)]);
        Gamma_i_per = Q_tmp(:,order+1:end)'; % Gamma orth
        M = Gamma_i_per*R_yf/R_uf;
        L = Gamma_i_per;
        left = zeros(ii*(no*ii-order), ni);
        right = zeros(ii*(no*ii-order), no*ii);
        for i = 1:ii
            left((i-1)*(no*ii-order)+1:i*(no*ii-order),:) = M(:,(i-1)*ni+1:i* ...
                                                              ni);
            right((i-1)*(no*ii-order)+1:i*(no*ii-order),1:(ii-i+1)*no) = ...
                L(:,(i-1)*no+1:no*ii);
        end
        right = right*[eye(no),zeros(no,order);zeros(no*(ii-1),no), ...
                       Gamma_im1];
        mat_db = right\left;
        mat_d = mat_db(1:no,:);
        mat_b = mat_db(no+1:no+order,:);
        if strcmp(method,'alg2')
            mat_a = pinv(Gamma_im1)*Gamma_ip1;
        elseif strcmp(method,'stable')
            mat_a = pinv(Gamma_i)*[Gamma_ip1;zeros(no,size(Gamma_i,2))];
        else
            error('method not supported!');
        end
    end
    
    system = idss(mat_a, mat_b, mat_c, mat_d); % write into matlab state space model
    

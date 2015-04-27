function hank_mat = build_hank(data, n_blk_r)
% build hankel matrix depend on the data matrix "data" and the block row number
% n_blk_r.
% data: the block vector, usually the rows are the data dimension, the colums
% are the sample points (may be time points)
% n_blk_r: number of block rows, the row dimension of the final hankel matrix
% is n_blk_r*data_r, where data_r is the number of rows of the block vector "data"
    
    [data_r, data_c] = size(data);
    n_c = data_c+1-n_blk_r; % determine the column number of the hankel matrix
    hank_mat = zeros(n_blk_r*data_r,n_c); % initiate                                   
    for i = 1:n_blk_r % place data
        hank_mat(data_r*(i-1)+1:data_r*i,:) = data(:,i:i+n_c-1);
    end
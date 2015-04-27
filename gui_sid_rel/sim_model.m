function [yh]=sim_model(T_start,T_step,t_0,t_end,app_end, offset,  yh, xh, At, Bt, Ct, Dt, p_vec, A1, B1, C1, D1)
A2t=reshape(At(:,:,2),size(A1)); B2t=reshape(Bt(:,:,2),size(B1)); C2t=reshape(Ct(:,:,2),size(C1)); D2t=reshape(Dt(:,:,2),size(D1)); 
A3t=reshape(At(:,:,3),size(A1)); B3t=reshape(Bt(:,:,3),size(B1)); C3t=reshape(Ct(:,:,3),size(C1)); D3t=reshape(Dt(:,:,3),size(D1)); 
A4t=reshape(At(:,:,4),size(A1)); B4t=reshape(Bt(:,:,4),size(B1)); C4t=reshape(Ct(:,:,4),size(C1)); D4t=reshape(Dt(:,:,4),size(D1)); 
A5t=reshape(At(:,:,5),size(A1)); B5t=reshape(Bt(:,:,5),size(B1)); C5t=reshape(Ct(:,:,5),size(C1)); D5t=reshape(Dt(:,:,5),size(D1)); 
A6t=reshape(At(:,:,6),size(A1)); B6t=reshape(Bt(:,:,6),size(B1)); C6t=reshape(Ct(:,:,6),size(C1)); D6t=reshape(Dt(:,:,6),size(D1)); 
A7t=reshape(At(:,:,7),size(A1)); B7t=reshape(Bt(:,:,7),size(B1)); C7t=reshape(Ct(:,:,7),size(C1)); D7t=reshape(Dt(:,:,7),size(D1)); 
A8t=reshape(At(:,:,8),size(A1)); B8t=reshape(Bt(:,:,8),size(B1)); C8t=reshape(Ct(:,:,8),size(C1)); D8t=reshape(Dt(:,:,8),size(D1)); 
A9t=reshape(At(:,:,9),size(A1)); B9t=reshape(Bt(:,:,9),size(B1)); C9t=reshape(Ct(:,:,9),size(C1)); D9t=reshape(Dt(:,:,9),size(D1)); 
A10t=reshape(At(:,:,10),size(A1)); B10t=reshape(Bt(:,:,10),size(B1)); C10t=reshape(Ct(:,:,10),size(C1)); D10t=reshape(Dt(:,:,10),size(D1)); 
A11t=reshape(At(:,:,11),size(A1)); B11t=reshape(Bt(:,:,11),size(B1)); C11t=reshape(Ct(:,:,11),size(C1)); D11t=reshape(Dt(:,:,11),size(D1)); 
for i=t_0:t_end
if(i==1) yh_m = sum(yh(:,1))/25+20;
else yh_m=sum(yh(:,i-1))/25+20; end
if(yh_m<T_start)
xh(:,i+1) = A1*xh(:,i) + B1*p_vec(:,i);
yh(:,i) = C1*xh(:,i) + D1*p_vec(:,i);
elseif(yh_m>T_start+T_step*(2-2) && yh_m< T_start+T_step*(2-1) )
xh(:,i+1) = A2t*xh(:,i) + B2t*p_vec(:,i);
yh(:,i) = C2t*xh(:,i) + D2t*p_vec(:,i);
elseif(yh_m>T_start+T_step*(3-2) && yh_m< T_start+T_step*(3-1) )
xh(:,i+1) = A3t*xh(:,i) + B3t*p_vec(:,i);
yh(:,i) = C3t*xh(:,i) + D3t*p_vec(:,i);
elseif(yh_m>T_start+T_step*(4-2) && yh_m< T_start+T_step*(4-1) )
xh(:,i+1) = A4t*xh(:,i) + B4t*p_vec(:,i);
yh(:,i) = C4t*xh(:,i) + D4t*p_vec(:,i);
elseif(yh_m>T_start+T_step*(5-2) && yh_m< T_start+T_step*(5-1) )
xh(:,i+1) = A5t*xh(:,i) + B5t*p_vec(:,i);
yh(:,i) = C5t*xh(:,i) + D5t*p_vec(:,i);
elseif(yh_m>T_start+T_step*(6-2) && yh_m< T_start+T_step*(6-1) )
xh(:,i+1) = A6t*xh(:,i) + B6t*p_vec(:,i);
yh(:,i) = C6t*xh(:,i) + D6t*p_vec(:,i);
elseif(yh_m>T_start+T_step*(7-2) && yh_m< T_start+T_step*(7-1) )
xh(:,i+1) = A7t*xh(:,i) + B7t*p_vec(:,i);
yh(:,i) = C7t*xh(:,i) + D7t*p_vec(:,i);
elseif(yh_m>T_start+T_step*(8-2) && yh_m< T_start+T_step*(8-1) )
xh(:,i+1) = A8t*xh(:,i) + B8t*p_vec(:,i);
yh(:,i) = C8t*xh(:,i) + D8t*p_vec(:,i);
elseif(yh_m>T_start+T_step*(9-2) && yh_m< T_start+T_step*(9-1) )
xh(:,i+1) = A9t*xh(:,i) + B9t*p_vec(:,i);
yh(:,i) = C9t*xh(:,i) + D9t*p_vec(:,i);
elseif(yh_m>T_start+T_step*(10-2) && yh_m< T_start+T_step*(10-1) )
xh(:,i+1) = A10t*xh(:,i) + B10t*p_vec(:,i);
yh(:,i) = C10t*xh(:,i) + D10t*p_vec(:,i);
elseif(yh_m> T_start+T_step*(11-2) )
xh(:,i+1) = A11t*xh(:,i) + B11t*p_vec(:,i);
yh(:,i) = C11t*xh(:,i) + D11t*p_vec(:,i);
end
end
yhs=yh;
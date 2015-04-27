function [yt] = Y_resp(n, t, t_offset, si, b, P)

if (n == 0)
    yt = 0;
    for i = 1 : length(si)
        yt = yt + (b(i)/si(i))*(exp(si(i)*(log(t)+t_offset))-1);
    end;
elseif (n == 1)
    yt = 10*Y_resp(0,t+1e10,t_offset,si,b,P) + P(1)*10*Y_resp(0,t,t_offset,si,b,P);
% elseif (n == 2)
%     if (P(n-1) == P(n))
%         yt = Y_resp(n-1,t+0.099,t_offset,si,b,P);
%     else
%         yt = Y_resp(n-1,t+0.099,t_offset,si,b,P) + P(n)*20*Y_resp(0,t,t_offset,si,b,P);
%     end
elseif (n > 1)
    if (P(n-1) == P(n))
        yt = Y_resp(n-1,t+0.1,t_offset,si,b,P);
    else
        yt = Y_resp(n-1,t+0.1,t_offset,si,b,P) + P(n)*20*Y_resp(0,t,t_offset,si,b,P);
    end
end
    
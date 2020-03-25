function Risk = RiskFunc(x, p)
% Calculate the failure risk of the ending time of each batch(comparing with the thresholds)
%
% INPUTS: 
%        x is predicted degradation path
%        p is production schedule information (positions for the ending time of each batch)
%
% OUTPUTS: 
%        Risks corespond to each order 
%
global Xf M N 
%% Calculate the failure risk of each datapoint
R_all = [];
for i = 1 : N
    n = 0;
    for j = 1:M
        if x(i,j)>= Xf
            n = n+1;
        else
            n = n;
        end
    end
    R_all = [R_all,n/M];
end
%% Calculate the failure risk of the ending time of production batch
% the position of the ending time
T = length(p);
R_batch = zeros(1,T); % T: number of batch
for i = 1: T
    R_batch(i) = R_all(p(i)); 
end
Risk = [];
Risk.all = R_all;
Risk.batch = R_batch;
end






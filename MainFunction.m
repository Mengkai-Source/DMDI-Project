%%%%%----- Simulated Annealing for discrete optimization (DMDII)----%%%%%

function [Xo, Xm, To, Tm, order_risk, Maintenance_Window] = MainFunction(MS, Health_Threshold, Cws2, ...
    Cwc2, Cm2, Cs2, L2, text, Cur_Date)
%
%  This is the main function to run the Dynamic Scheduler and display the preliminary results.
%
%  Inputs:
%       MS - represents HI_Curves indicated the prediction degradtion paths 
%       Health_Threshold - red threshold for the selected asset, used for the NU dynamic scheduler (highest threshold in deploy)
%       Cws2 - single value, spare part inventory cost for the selected asset
%       Cwc2 - single value, machine downtime cost for the selected asset
%       Cm2 - maintenance cost
%       Cs2 - shipping cost
%       L2 - single value, spare part lead time for the selected asset
%       text - cell matrix, production schedule info
%       Cur_Date -  current time indicating at which the prediction is
%       triggerd
%
%  Format of inputs:
%       text - An excel spreadsheet in '.xlsx' format contains information
%       in association to a specific machine (From left to right): 
%       Order -  Order ID (Column 1)
%       Operation - Type of operation(Column 2)
%       Material - Material needed for the operation(Column 3)
%       Scheduled Starting Time -  Scheduled Time when the order starts(Column 4)
%       Scheduled Finishing Time - Scheduled Time when the order finishes(Column 5)
%       Workplace -  Location where the operation is performed(Column 6)
%       Target Quantity - Amount of products of the order(Column 7)
%       Start Time -  When the order really starts(Column 8)
%       End Time - When the order really ends(Column 9)
%       Status - Machine in operation/not(Column 10)
%
%  Outputs:
%        Xo is the optimal threshold for spare part order placement
%        Xm is the optimal threshold for performing maintenance 
%        To is the countdown time for order placement 
%        Tm is the countdown time for maintenance performance
%        order_risk contains Risks coresponding to each order   
%        Maintenance_Window is the recommended maintenance opportunity
%        window suggesting early and late maintenance times
%
%  Format of outputs:
%        Current format displays two values of order placement threshold
%        and maintenance performance threshold and displays the corresponding
%        countdown time. In addition, actionable signal indication which can
%        be refered to suggest the actions with
%        respect to order purchase and maintenance performance.
%
%  NU_Version1.0
%  by Anqi He, Northeastern University, Boston

% Global
global X L Xf C1 C2 C3 Cws Cwc N M Cs C0 delta_t xlower xupper xdelt;

% Import the data of the future predicted degradation path
X = MS;

%% Define the feasible region

Xf = Health_Threshold;  % Provided by PDX-Deploy
L = L2; % Lead Time
xdelt= 0.005; %0.005
% initialize outputs
Xo = nan;
Xm = nan;
To = nan;
Tm = nan;
order_risk = nan;
% Pre-check the validity of the predicted paths
% Here, we need to guarantee the lower bound value of provided degradation
% paths all above the failure threshold
if length(find(X(end,:)<Xf))/size(X,2) < 0.8
    X(:,find(X(end,:)<Xf))=[];
    
    [row,col] = size(X);
    N = row; % The time length of predicted degradation paths( provided by PDX )
    M = col; % Number of predicted paths
    
    
    % Initial Parameters Setting
    C1 = Cm2; % Maintenance cost for type 1
    C2 = Cm2; % Maintenance cost for type 2
    C3 = Cm2; % Maintenance cost for type 3
    % Cws = 1; % Cost of inventory per day
    % Cwc = 1000; % Cost of machine breakdown per day
    Cws = Cws2;
    Cwc = Cwc2;
    Cs = Cs2;
    C0 = 0.5; % Initial cost can be removed
    delta_t = 1 ; % Unit time (provided by PDX)
    Pb = nan; p = nan;
    %% SA for Optimization
    % Create vectors containing the lower bound (|lb|) and upper bound constraints(|ub|).
    % Transform the bounds into discrete variables.
    xlower = min(min(X)); %%% X(1,1)min(min(X))
    xupper = Xf;
    l = [1 1];
    u = [length(xlower:xdelt:xupper) length(xlower:xdelt:xupper)];
    ini=ceil((u-l).*abs(rand(1,2)));%Initial guess for the minimum (threshold to order, threhold to perform maintenance)
    % Call the SA function to solve the problem with discrete variables.
    C=@CostRate;
    [x0,f0]=sim_anl(C,ini,l,u,400);
    
    % _Analyze the Results_
    xbestDisc = MapVariables(x0);
    
    Xo=xbestDisc(1); 
    Xm=xbestDisc(2); 
    
    %% Calculate the risk
    % Import the production schedule information
    % Extract the position of the ending time of each time points
    %%%%%%%%%%%%% Mengkai
    
    Pro_Schedule = text(2:end,4:5); % Extract the production schedule
    [r, c] = size(Pro_Schedule);
    for i = 1:c
        for j = 1:r
            Pro_Schedule{j,i} = datestr(Pro_Schedule{j,i}, 'mm/dd/yyyy HH:MM:SS'); %Convert the production date to 'mm/dd/yyyy HH:MM:SS'
        end
    end
    
    
    %Identify the current time corresponding to the production schedule
    if r>0
        for i =1:r
            Cur_base(i) = Pro_Schedule{i,2} - datetime(Cur_Date);
        end
        Cur_base.Format='d';
        positiveIndexes = Cur_base > 0;
        smallestPosValue = min(Cur_base(positiveIndexes));
    else
        smallestPosValue = [];
    end
    
    if ~isempty(smallestPosValue)
        rowOfSmallestValue = find(Cur_base == smallestPosValue);
        
        % Calculate the end time for each order
        Pro_Schedule = Pro_Schedule(rowOfSmallestValue:end,:);
        [ur, ~] = size(Pro_Schedule);
        T1 = datetime(Cur_Date);
        for i =1:ur
            T2 = datetime(Pro_Schedule{i,2});
            P(i) = T2 - T1;
            PB(i) = datetime(Pro_Schedule{i,1}) - T1;
        end
        PB.Format = 'd'; P.Format='d';  % Format the order start & end time duration in days
        Pb = round(days(PB));P=days(P);% The derived end ponits for each order
        p=round(P); % Convert duration array to numeric array by rounding the numerical value to the nearest integer
        
        % Probably should make NaNs, but that might break deploy, anyways the 
        % zeros at the end mean those schedules are more than 200 days out so they won't ever be visible anyways
        finalOrderRisk = [1:length(p);zeros(size(p))]; 
        if p(end)> row
            trunc_positive = row - p >0;
            smallest_trunc = min(trunc_positive(trunc_positive));
            ind = find(trunc_positive == smallest_trunc);
            p = p(1:max(ind));
            r = length(p);
        else
            ind = 1:length(p);
            p = p;
        end
       % Calculate risk for each order
        Risk = RiskFunc(X, p); 
        
        %% Dispaly risk corrsponding to each order
        order_risk=Risk.batch;
        
        order_ind =1:length(order_risk);
        order_risk=[order_ind;order_risk];
        
        finalOrderRisk(2,ind) = order_risk(2,:);
        order_risk = finalOrderRisk;
        
        %% Dispaly risk for each upcoming day
        Day_risk=Risk.all;
        Day_risk=Day_risk(2:end);
           
    end
    %Calculate To
    TO=[];
    for i=1:col
        ro=min(find(X(:,i)>=Xo));
        TO=[TO ro];
    end
    To=sum(TO)/length(TO)-1; % Countdown time from today 1 to the order palcement date To;
    To = floor(To);
    % countdown time array from the beginning time of prediction to order time
    
    %Calculate Tm
    TM=[];
    for i=1:col
        rr=min(find(X(:,i)>=Xm));
        TM=[TM rr];
    end
    Tm=sum(TM)/length(TM)-1; % Countdown time from today 1 to the maintenance date Tm;
    Tm = round(Tm);
    
    % Calculate the maintenance time window
    
    Pb = [Pb' p'];
    Indb = max(find(Pb(:,1) <= Tm));
    if Pb(Indb,2) >= Tm
        Maintenance_Window = Pb(Indb,:);
    elseif Indb < r
        Maintenance_Window = [Pb(Indb,2) Pb(Indb+1,1)];
    else
        Maintenance_Window=nan;
    end
else
    order_risk(1, 1:size(text,1) - 1) = 1:size(text,1) - 1;
    order_risk(2, 1:size(text,1) - 1) = 0;
end
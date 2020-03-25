%David Siegel
%April 28, 2018
%
% Modified By:      Matt Buzza
% Date Modified:    12/4/2018 (Uses the Unified_Dash_Scheduler function)
% Date Modified:    5/2/2019 Uses production schedule data now including
%                   updated MainFunction and Unified Dash Scheduler
% Purpose:          Simulates or loads data, and shows the output of the
%                   Unified_Dash_Scheduler and the predicted health curves 
%                   (the health curves don't need shown in the dashboard)
%  Outputs:
%          Xo is the optimal threshold for spare part order placement
%          Xm is the optimal threshold for performing maintenance 
%          To is the countdown time for order placement 
%          Tm is the countdown time for maintenance performance
%% Lets clear stuff
clear
clc
close all
%% Input data and parameters (These are from the unified dashboard, pull these from 
%the deploy dash,these should be the default settings)

% When do you want the trigger point to be? (first trigger point, or latest
% health) case1= 'firstTrigger', case2 = 'latestHealth'
triggerCase = 'latestHealth';

% What data do you want to load? originalData, simulatedData
dataCase = 'simulatedData2';

Health_Threshold = 1; % Highest threshold from deploy (red line)

% Model Inputs:
% prediction horizon (this should be adjustable for the application)
N_Ahead=200;

% Criteria used to trigger the health prediction: health index cutoff,
% correlation (trending) cutoff, and the window size for smoothing the 
% health index before applying the health index cutoff. We may want to have
% these options adjustable.
HI_Cut=0.35;
Corr_Cut=0.7;

% All default values
Cws = 1; % Inventory cost per unit time  %%% original 1
Cwc = 1000; % Machine downtime cost per unit time 10 %%% original 1000
L = 6;  % Lead time(settle down regarding the real case )
Cm = 25000; % Maintenance Cost
Cs = 500; % Shipping Cost

% Load Production Schedule Data (pull this prodSchedule from Forcam APIs)
[~,prodSchedule,~]= xlsread('C:\Users\xu.meng\Desktop\Dynamic Scheduler Test Code 20190908\Dynamic Scheduler Test Code 20190908\Data Files\Orders_Scheduled_01.01_MB.xlsx');
prodSchedule(2:end,:) = [];
dataFolderPath = fullfile(cd,'Data Files');

%% Simulate or Load Data (Ignore this part, just grab this same data from sandbox)

switch dataCase
    case 'originalData'
        % Load some example files...................................................
        %health value information
%         load('Z:\User\buzza\Lockheed Martin\MxD - Phase 1\Code\Dynamic Scheduler Test Code 20190702\Data Files\infoStream_5.mat',...
%             'healthValueAll');
        load(fullfile(dataFolderPath, 'infoStream_5.mat'),'healthValueAll')
        
        %Unit ID information
%         load('Z:\User\buzza\Lockheed Martin\MxD - Phase 1\Code\Dynamic Scheduler Test Code 20190702\Data Files\infoStream_3.mat',...
%             'unitID', 'timeStamp');
        load(fullfile(dataFolderPath, 'infoStream_3.mat'), 'unitID', 'timeStamp')
        
        % Pick Unit 2 as an example
        %lets user unit 10 as an example
        Unit_ID=10;
        HI_Temp = healthValueAll(unitID==Unit_ID);
        
        % Fix the timestamp so it aligns with the mock production data
        TS_Temp = 737243:1:737243+length(HI_Temp)-1;
        
    case 'simulatedData'
%         load('Z:\User\buzza\Lockheed Martin\MxD - Phase 1\Code\Dynamic Scheduler Test Code 20190702\Data Files\Simulated Spindle Health Values.mat')
        load(fullfile(dataFolderPath, 'Simulated Spindle Health Values.mat'))
        HI_Temp = HI;
        % Fix the timestamp so it aligns with the mock production data
        TS_Temp = datenum(2018,10,1):datenum(2018,10,1)+length(HI_Temp)-1;
        
    case 'simulatedData2'
%         load('Z:\User\buzza\Lockheed Martin\MxD - Phase 1\Code\Dynamic Scheduler Test Code 20190702\Data Files\20190801 Simulated Spindle Health Values.mat')
        load(fullfile(dataFolderPath, '20190801 Simulated Spindle Health Values.mat'))
        HI_Temp = healthValues;
        TS_Temp = timestamp';
        
        % Just include the data up until 5/30/2019, that is when the Tm
        % suddenly dropped to zero
        iKeep = find(TS_Temp < datenum(2019, 5, 31), 1, 'last');
        HI_Temp = HI_Temp(1:iKeep);
        TS_Temp = TS_Temp(1:iKeep);
        
    otherwise
        error('Please input an acceptable dataCase')
end

%% Depending on user input, only input the data up until the first trigger 
% point, or input all the data (if you input all the data, the
% Unified_Dash_Scheduler will just check if the last health value meets the
% trigger criteria, and then run the scheduler to calculate To, Tm, and the
% order risks)
switch triggerCase
    
    case 'firstTrigger'
        % Trigger criteria just looks at the last Nw health values
        Nw = 50;
        
        % Initialize
        Trigger = zeros(length(HI_Temp),1);
        
        % Store the trigger values, for testing purposes
        stored_values = zeros(length(HI_Temp), 2);
        
        %find any trigger points
        for ii=Nw:1:length(HI_Temp)
            This_HI_Temp = HI_Temp(ii-Nw+1:ii,1);
            This_TS_Temp = TS_Temp(ii-Nw+1:ii)';
            
            stored_values(ii, 1) = corr(This_TS_Temp,This_HI_Temp);
            stored_values(ii,2) = This_HI_Temp(end);
            
            if ii >= Nw && This_HI_Temp(end)>HI_Cut && stored_values(ii, 1)>Corr_Cut
                Trigger(ii,1)=1;
            else
                Trigger(ii,1)=0;
            end
        end
        
        HI_Temp = HI_Temp(1:find(Trigger,1,'first'));
        TS_Temp = TS_Temp(1:find(Trigger,1,'first'));
        
    case 'latestHealth'
        % Don't do anything, just use all the HI data
    otherwise
        error('Please Input an acceptable triggerCase')
end

%% Run scheduler (This contains PDX prediction code and NUs dynamic 
% scheduler code that calculates Xo,Xm,To,Tm,order_risk which are shown in
% the dashboard) RUN THIS EVERYTIME THERE IS A NEW HEALTH VALUE CALCULATED

% Make sure later to remove orders that ended before the latest health
% value, we only want to show current order and future orders in the 
% dashboard anyways
prodTimes = nan(length(prodSchedule(:,4))-1, 2);
for rowInd = 2:length(prodSchedule(:,4))
    prodTimes(rowInd-1, 1) = datenum(prodSchedule{rowInd,4}, 'mm/dd/yyyy HH:MM:SS PM');
    prodTimes(rowInd-1, 2) = datenum(prodSchedule{rowInd,5}, 'mm/dd/yyyy HH:MM:SS PM');
end
rInd = find(prodTimes(:,2) < TS_Temp(end));
prodTimes(rInd,:) = [];
prodSchedule(rInd+1,:) = [];


% Calculates Xo,Xm,To,Tm,order_risk and maintenance opportunity window
if size(TS_Temp,1) < size(TS_Temp,2)
    TS_Temp = TS_Temp';
end
[Xo, Xm, To, Tm]...
    = Unified_Dash_Scheduler_v3(HI_Temp,TS_Temp,Cws, Cwc, Cm, Cs, L, N_Ahead,...
    HI_Cut, Corr_Cut, Health_Threshold, prodSchedule);

%% Plot the generated HI_Curves (Ignore this part, just for some validation)

% Remove NaNs from HI_Temp (Normally done inside Unified_Dash
R_Ind = isnan(HI_Temp);
HI_Temp(R_Ind) = [];
TS_Temp(R_Ind) = [];

% Northeasterns algorithm will ignore production orders that happen more
% than N_Ahead samples in the future
prodTimes = prodTimes(1:size(order_risk,2),:);

ylims = [0 1.1];

% Green to red color gradient
colorsGR = jet(20);
colorsGR(1:10,:) = [];

% Only plot if the scheduler was triggered
if ~isnan(Xo)
    
    [HI_Curves, Trigger_Point] = ...
        Prediction_Routine_v9(HI_Temp, TS_Temp, N_Ahead, HI_Cut, Corr_Cut);
    
    
    % Done inside Northeasterns function
    % Pre-check the validity of the predicted paths
    % Here, we need to guarantee the lower bound value of provided degradation
    % paths all above the failure threshold
    HI_Curves(:,HI_Curves(end,:)<Health_Threshold)=[];
    HI_Curves(:,HI_Curves(1,:)>Health_Threshold) = []; % Remove curves that are already above the threshold
    
    % Plot the prediction curves (for demonstration)
    %make a plot and go from their
    Index=1:1:Trigger_Point;
    
    Warning=0.5;
    Failure=1;
    
    figure
    % Plot all the health values from the last baseline
    plot(TS_Temp(1:1:Trigger_Point),HI_Temp(1:1:Trigger_Point),'b-o')
    grid on
    hold on
    TS = [TS_Temp(1:Trigger_Point); [TS_Temp(Trigger_Point)+1:TS_Temp(Trigger_Point)+N_Ahead]'];
    
    % Plot the thresholds from deploy
    plot([TS(Index(1)) TS(Index(end)+N_Ahead)],[Warning Warning],'y','linewidth',2)
    plot([TS(Index(1)) TS(Index(end)+N_Ahead)],[Failure Failure],'r','linewidth',2)
    
    % Plot the production schedule and the risk for each batch
    for rowInd = 1:size(prodTimes,1)
        x = [prodTimes(rowInd,1) prodTimes(rowInd,2) prodTimes(rowInd,2)...
            prodTimes(rowInd,1)];
        y = [ylims(1) ylims(1) ylims(2) ylims(2)];
        
        % Get the color of the batch, based on order risk
        ind = ceil(order_risk(2, rowInd)*10);
        if ind == 0
            ind = 1;
        end
        p=patch(x,y,colorsGR(ind,:));
        set(p,'FaceAlpha',0.5);
        text(prodTimes(rowInd,1), .1, num2str(round(order_risk(2, rowInd),3)))
    end
    
    datetick('x', 'mm/dd/yy')
    xlim([TS(Index(1)) TS(Index(end)+60)])
    ylim(ylims)
    
for CurveInd = 1:size(HI_Curves,2)
    x = TS(Trigger_Point:N_Ahead+Trigger_Point-1);
    y = HI_Curves(1:N_Ahead, CurveInd);
    y(end) = NaN;
    patch(x,y,[.85 .85 .85], 'EdgeAlpha',.01,'linewidth',4)
end

    xlabel('Cycle #','fontsize',14)
    ylabel('HI','fontsize',14)
    title(['Health Index over Time To=' num2str(To) ' Tm=' num2str(Tm)],...
        'fontsize',14)
    legend('HI','Warning','Failure','location','best');
    
end
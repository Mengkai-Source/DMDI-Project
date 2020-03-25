function [Xo, Xm, To, Tm, order_risk,Maintenance_Window] = ...
    Unified_Dash_Scheduler_v3(HI_Temp, TS_Temp, Cws, Cwc, Cm, Cs, L, N_Ahead,...
    HI_Cut, Corr_Cut, Health_Threshold, prodSchedule)
% Author:       Matt Buzza
% Date Created: 11/26/2018
% Purpose:      Takes all the health values since the last baseline,
%               triggers NUs scheduler if the health is low and has a good
%               trend. Run this if new health values have been calculated,
%               or Cws, Cwc, L, N_Ahead, HI_Cut, Corr_Cut, Health_Threshold
%               parameters have been changed in the UI
%
% Date Modified: 1/21/2019 v2 - Added timestamp as an input to this
%                               function. Also updated it so it first takes
%                               the average health value from each day for
%                               smoothing. The rest of the code won't work
%                               properly if there is more than 1 health
%                               value per day. Also changed it so it only
%                               runs Northeasterns Scheduler part if To and
%                               Tm are not already calculated. It also does
%                               the countdown for To and Tm, each day there
%                               is a new data record it will countdown To
%                               and Tm by 1.
% Date Modified: 5/2/2019 v3 - Incorporates the production schedule, new
%                              input. Now Uses NUs MainFunction to
%                              calculate To, Tm, and order risk everytime
%                              new health values are available.
%
%
%   Inputs:
%   HI_Temp - nx1 double array, all the health values since the last 
%             baseline for the selected asset
%   TS_Temp - Date (format: num) represents the historically temporal operation on
%   calendar
%   Cws - single value, spare part inventory cost for the selected asset
%   Cwc - single value, machine downtime cost for the selected asset
%   Cm - maintenance cost
%   Cs - shipping cost
%   L - single value, spare part lead time for the selected asset
%   N_Ahead - "Prediction horizon used for scheduler"
%   HI_Cut - "Health index cutoff used for triggering scheduler"
%   Corr_Cut - "Correlation cutoff used for triggering scheduler"
%   Health_Threshold - red threshold for the selected asset, used for the
%                      NU dynamic scheduler (highest threshold in deploy)
%   prodSchedule - cell matrix, production schedule info
%
%
%   Outputs:
%   Xo - threshold for ordering the parts set by Northeasterns  scheduler (MainFunction.m)
%   Xm - threshold for when to do maintenance set by Northeasterns  scheduler (MainFunction.m)
%   To - predicted time remaining for when Xo will be reached set by Northeasterns  scheduler (MainFunction.m)
%   Tm - predicted time remaining for when Xm set by Northeasterns  scheduler (MainFunction.m)
%   order_risk - risk of failure before the end of each order set by Northeasterns  scheduler (MainFunction.m)
%   Maintenance_Window is the recommended maintenance opportunity window suggesting early and late maintenance times(MainFunction.m)


% Remove NaNs from HI_Temp
R_Ind = isnan(HI_Temp);
HI_Temp(R_Ind) = [];
TS_Temp(R_Ind) = [];


[HI_Curves, Trigger_Point] = Prediction_Routine_v9(...
    HI_Temp, TS_Temp, N_Ahead, HI_Cut, Corr_Cut);

    
% Check if scheduler was triggered (will be NaN if it wasn't)
if isnan(HI_Curves)
    % Set the outputs to NaN if output wasn't triggered
    Xo = nan;
    Xm = nan;
    To = nan;
    Tm = nan;
    
    % Check to see if there are even any current or future orders
    if size(prodSchedule,1) > 1
        % Set all the order risks to 0 if scheduler wasn't triggered
        order_risk(1, 1:size(prodSchedule,1) - 1) = 1:size(prodSchedule,1) - 1;
        order_risk(2, 1:size(prodSchedule,1) - 1) = 0;
    else
        order_risk = nan;
    end
    
    % Check if the health is already above the Maintenance threshold, if it is
    % then make Xo=Xm=Threshold, and To=Tm=0, and order risk all 1's
elseif HI_Temp(Trigger_Point) > Health_Threshold
    Xo = Health_Threshold;
    Xm = Health_Threshold;
    To = 0;
    Tm = 0;
    
    % Check to see if there are even any current or future orders
    if size(prodSchedule,1) > 1
        order_risk(1,1:size(prodSchedule,1) - 1) = 1:size(prodSchedule,1) - 1;
        order_risk(2,1:size(prodSchedule,1) - 1) = 1;
    else
        order_risk = nan;
    end
    
    % If scheduler was triggered, the health isn't already above the
    % threshold, and To and Tm are not already calculated then run the 
    % scheduler
else
    % Remove HI_Curves that start above the threshold
    R_Ind = HI_Curves(1,:) > Health_Threshold;
    HI_Curves(:,R_Ind) = [];
    
    
    % Cur_Date is the date of the trigger point, which is also the date of
    % the last health value
    Cur_Date = datestr(TS_Temp(Trigger_Point), 'mm/dd/yyyy HH:MM:SS');
    
    [Xo, Xm, To, Tm, order_risk,Maintenance_Window] = MainFunction(HI_Curves, Health_Threshold, Cws, ...
        Cwc, Cm, Cs, L, prodSchedule, Cur_Date);
    
end
end


function [HI_Curves, Trigger_Point] = ...
    Prediction_Routine_v9(HI_Temp, TS_Temp, N_Ahead, HI_Cut, Corr_Cut)

% Modified By:      Matt Buzza
% Date Modified:    11/1/2018 (Make some of the settings function inputs, 
%                   rather than hard coded, in case we want to make them
%                   adjustable from the unified dashboard.)
% Date Modified:    11/26/2018 (Remove NaNs from HI_Temp, to avoid errors
%                   if this ever occurs. Also, returns HI_Curves as NaN if
%                   there was no triggering. Hard code Nw now. Output
%                   HI_Smooth. HI curves are generated from smoothed data.
%                   Smoothed data is done using median fitler now.)
% Date Modified:    1/21/2019 v8 - Added time stamp as an input and output.
%                   It now takes the average health values for each day, 
%                   before doing anything else. Since previously the code 
%                   won't work properly if there are more than one health 
%                   value per day. 
% Date Modified:    5/6/2019 v9 - Switched to a linear fit. Uses last 50 HVs
%                   instead of every HV. Changed it so timestamp is the
%                   predictor variable, since there are many sources of
%                   gaps in the health values (IPC shutoff, MTConnect
%                   messing up, operator forgetting to run the routine
%                   program). Changed it so the trigger point is always the
%                   last health value, since that is how it will be when
%                   the system is deployed. Made it do robust regression,
%                   so smoothing is not necessary. No longer need to take
%                   mean health value for each day, since timestamp is now
%                   predictor variable. Set a random seed. Only checks
%                   correlation of the last Nw health values
%
%
%
%   Inputs:
%   HI_Temp - mx1 array of health indices
%   TS_Temp - mx1 array of time stamps
%   N_Ahead - single value, prediction horizon
%   HI_Cut - single value, health index cutoff
%   Corr_Cut - single value, correlation cutoff
%
%   Outputs:
%   HI_Curves - mmxn matrix of prediction curves
%   Trigger_Point - if the prediction was triggered for last HV, the trigger point is
%                   is just the index of the last health value. Otherwise
%                   the trigger point is 0, and the dynamic scheduler
%                   (MainFunction) wont be called


%Trigger Logic (make it relatively simple to start)

%Criteria 1 (HI>HI_Cut) %at least 50% degradation

%Criteria 2 (Correlation with health index >Corr_Cut) - noticable trend

%Criteria 3 (At least 50 previous samples) - enough data to do the
%prediction

% Focus on last 50 health values
Nw = 25;

%%

% Make sure the last point was a trigger point, otherwise set the trigger
% as empty/not triggered
if length(HI_Temp) >= Nw
    corrCheck = corr(TS_Temp(end-Nw+1:end), HI_Temp(end-Nw+1:end));
    if HI_Temp(end) >= HI_Cut && corrCheck >= Corr_Cut
        Trigger_Point = length(HI_Temp);
    else
        Trigger_Point = [];
    end
else
    Trigger_Point = [];
end



% If there is no trigger, set HI_Curves as NaN
if isempty(Trigger_Point)
    HI_Curves = nan;
    
else

    % Just use the last Nw for prediction
    Index=Trigger_Point-Nw+1:1:Trigger_Point;
    HI_Trend = HI_Temp(Index);
    TS_Trend = TS_Temp(Index);
    
    % linear fit
    fitresult = fit(TS_Trend,HI_Trend,'poly2','Robust','Bisquare');
    p1 = fitresult.p1;
    p2 = fitresult.p2;
    p3 = fitresult.p3;
    
    %Index2
    Index2 = (1:N_Ahead) + TS_Temp(end); % Predict in terms of days
    
    Ypred = p1*Index2.^2 + p2*Index2 + p3';
    Ypred = Ypred';
    
    % prediction intervals
    alpha=0.05; %95 confidence level
    ci = predint(fitresult,Index2,1-alpha);
    
    % Curves = nan * ones(N_);
    STDEV_Prediction=ci(:,2)-Ypred;
    N_Distribution=1000; %1000 samples
    rng(1)
    for bbb=1:N_Distribution
        Curves(bbb,:)=randn(1,1)*STDEV_Prediction+Ypred;
    end
    
    
    %health value distribution
    HI_Curves=Curves';
    
end


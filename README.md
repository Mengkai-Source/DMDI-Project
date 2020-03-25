Dynamic-Scheduler-with-Production-Schedule
This version with SA and incorporate the production schedule

Dynamic-Scheduling-SA Prof.Xiaoning Jin, Prof.Kamarthi, PhD student: Anqi He, Mengkai Xu

This is a preliminary demo of the dynamic scheduler for test implementation. There are six Matlab function files (MATLAB code) in the folder. One sample data file : 'path.mat', and five Matlab functions: 'MainFunction.m', 'CostRate.m', 'sim_anl.m', 'mu_inv.m', 'MapVariables.m', 'RiskFunc.m'. Per the discussion with PDX and FORCAM, Current format displays the countdown time To (Order placement) and Tm (Maintenance performance), a vector contains the risk for to the current time and orward, a vector contains the risk corresponding to each order in the production shcedule.

#1. Input & Output

Input Information Input Data file: Prediction Dataset: Sample_Data (Cur_Date, HI_Curves, Mach_ID); Production Schedule Dataset:P Cur_Date: Date when prediction is triggered HI_Curves: N X M matrix which consists of multiple degradation paths of predicted machine health values from PDX N: Time index M: Path index

The machine health values of the first row should be same because all paths start at the same prediction trigger time Mach_ID: Machine ID (workplace in the production shcedule spreadsheet) Production Schedule Dataset: machine production shcedule, i.e. start working time, and end working time.

Output Information Output Format Risk.all: numeric array in which each number represents the failure risk of each predicted time. Risk.batch: numeric array in which each number represents the failure risk of each scheduled production batch (these values need to compare with some thresholds to determine the warning color of each batch) Countdown time for order: numeric array Countdown time for maintenance: numeric array

#2. Runing environment: MATLAB (No version requirement & No Toolbox is needed)

#3. Function file: Install

'MainFunction.m':This is the main function to run the Dynamic Scheduler and display the desired results. The Outputs from running this function contain: Countdown time To (Order placement) and Tm (Maintenance performance); A vector contains the risks corresponding to current time onwards; A vector contains the risks corresponding to each order in the production schedule. (Thress adjustable parameters: Cws (Spare part inventory cost Dollar ($)/ per day), Cwc (Machine downtime cost Dollar ($)/ per day) and L (Part lead time: Days). The parameter for lead time is 'L'.)

'CostRate.m': This function is the objective function to calculate the 'cost' metric based on the predicted degradation path.

'MapVariables.m': MapVariables interatively assign values to the elements in vector x

'sim_anl.m': Minimizes a function with the method of simulated annealing

'mu_inv.m': This function is used to generate new point according to lower and upper and a random factor proportional to current point.

'RiskFunc.m': This function is used to calculate the risks.

#4. Run

Run the main program: "MainFunction"

Note: Please make sure the date for the loaded prediction data matches the date for the production schedule.

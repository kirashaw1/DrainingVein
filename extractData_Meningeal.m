% Script to extract the time series metrics across individual walking trials 
% for the draining vein awake mouse data into tables for stats in R. 
% Written by Kira, June 2024

% functions needed: findDeltaF, getTraceParam 

clear all; close all; 

%% SPECIFY PREFS for functions

% Time series metrics prefs (for getTraceParam)
%specify which time period you want to look for your ts parameters
%between, i.e. look for AUC and max peak
%reqTim means requested time
%if this is left empty, some of the parameters cannot be detected and
%will be sent out as empty variables
prefs.reqTim = [5 10]; %[] %seconds
%the time (secs) you want to average either side of your detected max
%peak, i.e. to make sure it isn't just a noise spike - dnt make this a
%large number (best to be <0.5)
%this will be used for all peaks detected
prefs.pkCheck = 0.25; %seconds
%if you want to look at the value of the trace at a certain time point
%(secs)
%reqPk means requested peak
%if this is left empty, some of the parameters cannot be detected and
%will be sent out as an empty variable
prefs.reqPk = []; %[] %seconds

%% data load and clean 
fps = 7.9932; %2DOIS sampling rate
time2norm = [1,40]; %baseline time to norm in frames (1st 5 seconds)

%load data: draining vein
load('locoHaemDat_perLocoEvent_meningeal57.mat')
%rename and normalise traces 
%NB exports trace as delta D/D, multiply by 100 if want as %
for a = 1:size(hbo_walk_tot,1)
    hbo_mn(a,:) = findDeltaF(hbo_walk_tot(a,:),time2norm);
    hbr_mn(a,:) = findDeltaF(hbr_walk_tot(a,:),time2norm);
    hbt_mn(a,:) = findDeltaF(hbt_walk_tot(a,:),time2norm);
end
clear hbo_walk_tot hbr_walk_tot hbt_walk_tot a; %clear old variables

%transpose animal ID to be same dimension as other variables
animalID = animalID'; 

% create time vector
time_haem = [0:size(hbo_mn,2)-1]/fps; %create time vector 

% resize walking trace to same length using interp1 function
%NB loco recording time is out from 2DOIS recording time by 0.000# seconds-
%is this okay? 
fps_loco = size(walk_tot,2)/time_haem(end); 
time_loco = [0:size(walk_tot,2)-1]/fps_loco;
for a=1:size(walk_tot)
    walk_tot_interp(a,:) = interp1(time_loco, walk_tot(a,:), ...
        linspace(0, time_loco(end), size(hbt_mn,2))); 
end
clear a

%% time series metrics 

% time series parameters- call function getTraceParam
for a = 1:size(hbt_mn,1) %loop loco event trials
    %draining vein
    [traceparameters_mn(a)] = getTraceParam(hbt_mn(a,:),fps,prefs);
    %loco
    [traceparameters_loco(a)] = getTraceParam(walk_tot_interp(a,:),fps,prefs);
end
clear a

%save matfile with traces
save([cd,filesep,'dvpaper_tstraces_meningeal.mat'],"time_haem","traceparameters_mn",...
    "traceparameters_loco","walk_tot_interp","fps","animalID","grpID","sessID","trialID",...
    "hbt_mn", "hbo_mn", "hbr_mn");

%% Create data table for stats
%transform data so each vessel category is row upon row for R table
%analysis
vessID=[]; maxpk=[]; t2max=[]; t2o=[]; minpk=[]; t2min=[]; loco_auc=[]; 
for b = 1:size(animalID,1)
    vessID{size(vessID,1)+1,1} = 'mn';
    maxpk(size(maxpk,1)+1,1) = traceparameters_mn(b).maxPk*100; %convert to percent
    t2max(size(t2max,1)+1,1) = traceparameters_mn(b).t2p;
    t2o(size(t2o,1)+1,1) = traceparameters_mn(b).t2o;
    minpk(size(minpk,1)+1,1) = traceparameters_mn(b).minPk*100; %convert to percent
    t2min(size(t2min,1)+1,1) = traceparameters_mn(b).t2min;
    loco_auc(size(loco_auc,1)+1,1) = traceparameters_loco(b).AUC;
end
clear b 

% %save table with metrics for stats in R 
% T = table(grpID, animalID, sessID, trialID, vessID, maxpk, t2max, t2o, ...
%     minpk, t2min, loco_auc);
% writetable(T,[cd,filesep,'data4LMM_meningeal.xlsx']);

figure;
hold on; 
plot(time_haem,mean(hbt_mn,1),'g');
plot(time_haem,mean(hbr_mn,1),'b');
plot(time_haem,mean(hbo_mn,1),'r');






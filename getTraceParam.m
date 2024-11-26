
function [traceparameters] = getTraceParam(trace,fps,prefs)

%function to send in a trace, and get some of the parameters out, e.g.
%maximum peak size, time of max peak, area under the curve, value of trace
%at requested time point (e.g. at stim off point), time to peak, time on
%onset (25% of peak)
% REQUIRES FUNCTION INPAINT_NANS - this is because trapz (for AUC) cannot
% work if there is a nan in the trace (just outputs NaN)- so can fill the
% nans in instead 

%INPUTS:
% trace - this should be the data trace, and should be sent in as a 2D
%         trace. Make sure you send it in as trials in dim1, and time in dim2.
%         Loop further dimensions outside of the function. NB/ better
%         to send in the full trace, i.e. not just the time period between
%         the peaks.
% fps   - the frames per second info (NB/ for a line scan send in the sps)
% prefs - these specify what timings you want for the parameters - see the
%         nargin < 3 info below, specify these outside the function.


%OUTPUTS:
% traceparameters - this is a structure, containing all the parameters
%                   inside. Access them by asking for traceparameters.auc
%                   for example.

if nargin < 3
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
    prefs.reqPk = [10]; %[] %seconds
end

% check if the trace variable is 2D
if ndims(trace) == 2
    
    % make a time vector to get info out in seconds
    if fps > 1 %xy movies
        time = [0:size(trace,2)-1]/fps;
    else %linescan
        time = [0:size(trace,2)-1]*fps;
    end %end of check if xy or linescan recording
    
%     % check data trace is in correct order
%     if size(trace,1) > size(trace,2) %check dimensions of trace sent in
%         %display a warning if 2nd dim is larger, in case accidentally sent in
%         %wrong way round (as time is usually a larger dim than trials, but not
%         %always)
%         disp('Warning, dim1 of trace larger than dim2. Make sure time is in dim2');
%     end %end of check trace dimensions order
    
    %% STEP 1: get parameters (AUC, t2p, pkSz, etc) within requested period
    %(e.g. during stim) - prefs.reqTim
    
    %extract time period requested for finding peaks etc
    time_ttt = time(:,find(time>=prefs.reqTim(1),1): ...
        find(time>=prefs.reqTim(2),1));
    %preference for checking either side of peak (in seconds) - convert to
    %frames
    %Note: you are limited by your Hz, so will not necessarily be exactly
    %accurate. The faster the Hz, the more accurate this will be.
    framesPkCheck = find(time>=prefs.pkCheck,1);
    
    for a = 1:size(trace,1) %loop trials from 2D data trace
        
        %extract trace during requested time period - create temp variable
        trace_ttt = trace(a,find(time>=prefs.reqTim(1),1): ...
            find(time>=prefs.reqTim(2),1));
        
        %% MAX PK PARAMETERS:
        
        % SIZE OF PK:
        %get out the maximum peak, and its index
        [~,pkInd_ttt] = max(trace_ttt);
        
        if ~isempty(pkInd_ttt)
            %convert this back to an index across the whole trace
            pkInd = pkInd_ttt + (find(time>=prefs.reqTim(1),1)-1);
            
            %check if there is enough time around the peak to add and subtract
            %these frames either side
            if pkInd > framesPkCheck && size(trace(a,:),2) >= pkInd+framesPkCheck
                %the max pk to save out is an average either side of the max peak
                %by the specified time period
                traceparameters.maxPk(a) = mean(trace ...
                    (a,pkInd-framesPkCheck:pkInd+framesPkCheck),2);
            else %not enough time
                %max peak is a NaN
                traceparameters.maxPk(a) = NaN;
            end %end of check if enough time
            
        else
            traceparameters.maxPk(a) = NaN;
        end
        
        %TIME TO PK:
        %this will be the time after the requested period (e.g. 5 seconds),
        %rather than the time from the very start of the trace (which includes
        %the baseline)
        %check if you have a max peak value (if a nan no pt getting t2p)
        if ~isnan(traceparameters.maxPk(a))
            if fps > 1 %xy movies
                traceparameters.t2p(a) = pkInd_ttt/fps;
            else %linescan
                traceparameters.t2p(a) = pkInd_ttt*fps;
            end %end of check if xy or ls recordings
        else %no peak val (NaN)
            traceparameters.t2p(a) = NaN;
        end %end of check if have a peak val
        
        %TIME TO ONSET:
        %this will be the first time pt which reaches >=25% of the max peak
        %value (i.e. to show the dilation has started kicking in)
        %NOTE this doesn't always work when there isn't a clearly increasing
        %dilation - if it is a noisy trace for example, it might already be at
        %25% in the first frame - so may not be worth presenting always?
        %depends on data
        if ~isnan(traceparameters.maxPk(a))
            if fps > 1 %xy movies
                if ~isempty(find(trace_ttt>=(trace_ttt(:,pkInd_ttt)/4),1))
                    traceparameters.t2o(a) = find(trace_ttt>= ...
                        (trace_ttt(:,pkInd_ttt)/4),1) / fps;
                else
                    traceparameters.t2o(a) = NaN;
                end
            else %linescan
                traceparameters.t2o(a) = find(trace_ttt>= ...
                    (trace_ttt(:,pkInd_ttt)/4),1) * fps;
            end
        else %no max peak
            traceparameters.t2o(a) = NaN;
        end %end of check if any max peak detected


        %get out the maximum peak, and its index
        clear pkInd_ttt; 
        [~,pkInd_ttt] = min(trace_ttt);
        
        if ~isempty(pkInd_ttt)
            %convert this back to an index across the whole trace
            pkInd = pkInd_ttt + (find(time>=prefs.reqTim(1),1)-1);
            
            %check if there is enough time around the peak to add and subtract
            %these frames either side
            if pkInd > framesPkCheck && size(trace(a,:),2) >= pkInd+framesPkCheck
                %the max pk to save out is an average either side of the max peak
                %by the specified time period
                traceparameters.minPk(a) = mean(trace ...
                    (a,pkInd-framesPkCheck:pkInd+framesPkCheck),2);
            else %not enough time
                %max peak is a NaN
                traceparameters.minPk(a) = NaN;
            end %end of check if enough time
            
        else
            traceparameters.minPk(a) = NaN;
        end

        %TIME TO MIN:
        %this will be the time after the requested period (e.g. 5 seconds),
        %rather than the time from the very start of the trace (which includes
        %the baseline)
        %check if you have a min peak value (if a nan no pt getting t2p)
        if ~isnan(traceparameters.minPk(a))
            if fps > 1 %xy movies
                traceparameters.t2m(a) = pkInd_ttt/fps;
            else %linescan
                traceparameters.t2m(a) = pkInd_ttt*fps;
            end %end of check if xy or ls recordings
        else %no peak val (NaN)
            traceparameters.t2m(a) = NaN;
        end %end of check if have a peak val
        
        %% AREA UNDER THE CURVE:
        
        % this will also be detected during the specified time period (e.g.
        % during the stim 5-10secs, or a bit longer, e.g. 5-12 secs) - it is
        % the same time period within which the time to peak is detected
        %check not too many NaNs in data- so infilling them is just "made
        %up data" per se
        if size(find(isnan(trace_ttt)),2) > size(trace_ttt,2)/2 %i.e. 50% nans
            traceparameters.AUC(a) = NaN; 
        else
            if isempty(find(isnan(trace_ttt)))
                %no nans
                traceparameters.AUC(a) = trapz(time_ttt,trace_ttt);
            else
                %nans to fill in 
                traceparameters.AUC(a) = trapz(time_ttt,inpaint_nans(trace_ttt));
            end
        end
        
        %% REQUESTED PEAK
        
        % this looks for the value of the data trace at the time pt (in
        % seconds) the user has specified. NB/ if the user doesn't need to look
        % at a certain time point, just leave this as an empty bracket and the
        % function won't look for this. This is sometimes useful when it looks
        % like there might be differences when the stim goes off for EG - so at
        % 10seconds.
        
        if ~isempty(prefs.reqPk) %check if requested pk check at certain time pt
            %NB/ still take the requested time (e.g. 0.25secs) either side of
            %requested frame, to erradicate noise spikes
            traceparameters.reqPk(a) = mean(trace(a, ...
                find(time>=prefs.reqPk,1)-framesPkCheck : ...
                find(time>=prefs.reqPk,1)+framesPkCheck));
        else %didn't request specific time pt
            traceparameters.reqPk(a) = NaN;
        end %end of check if specified a peak check at a certain time point
        
        %clear the variables which you will use again in the loop
        clear trace_ttt pkInd pkInd_ttt;
        
    end %end of looping trials from 2D data trace
    
else %sending in a 3d var - need to loop outside trace
    
    disp('too many dims in trace, loop outside func. Exiting...');
    
end %end of warning message


end %end of function


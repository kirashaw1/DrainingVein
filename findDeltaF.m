function [dataNorm]=findDeltaF(data2norm, frames2norm)

%function written by Kira, August 2017
%function to normalise data by a specified baseline period - uses delta F/F
%note that the time dimension (to normalise) is assumed as the last dim

%INPUTS-
%data2norm   : send in the data to norm (can be 2 or 3 dims, time dim last)
%frames2norm : frames within which to normalise data (i.e. not in secs)
%OUTPUTS-
%dataNorm    : sends out the data now normalised to zero (within frames
%              specified)


%find out if the data has 2 or 3 dimensions
t=ndims(data2norm);

%determine if data is a cell (i.e. when averaging across ALL sessions)
if iscell(data2norm)==0 %inputted data is not a cell
    %for 2D data, just need to loop through 1st dim, and norm
    if t==2
        for i=1:size(data2norm,1) %loop dim 1
            %subtract baseline, then divide by it...
            %abs the divided base (but not subtracted)
            %tested this, only way preserves shape around zero and not neg
            if size(frames2norm,2)==2 %sent in just one norm period
               dataNorm(i,:) = (data2norm(i,:)-(nanmean(data2norm(i,...
                    frames2norm(1):frames2norm(2)))))./abs(nanmean(data2norm...
                    (i,frames2norm(1):frames2norm(2))));
            else %sent in two normalisation periods
                dataNorm(i,:) = (data2norm(i,:)-(nanmean(data2norm(i,...
                    [frames2norm(1):frames2norm(2), ...
                    frames2norm(3):frames2norm(4)]))))./abs(nanmean(data2norm...
                    (i,[frames2norm(1):frames2norm(2),...
                    frames2norm(3):frames2norm(4)])));
            end %end of check if have multiple norm periods
        end %end dim 1 loop
        %for 3D data, need to loop through 1st 2 dims and norm
    elseif t==3
        %normalise data by dividing by the mean of the normalisation period
        for i=1:size(data2norm,1) %loop dim 1
            for j=1:size(data2norm,2) %loop dim 2
                if size(frames2norm,2)==2 %sent in just one norm period
                    dataNorm(i,j,:)= (data2norm(i,j,:)-(nanmean(squeeze(...
                        data2norm(i,j,frames2norm(1):frames2norm(2))))))...
                        ./abs(nanmean(squeeze(data2norm(i,j,frames2norm(1):...
                        frames2norm(2)))));
                else %sent in two normalisation periods
                   dataNorm(i,j,:)= (data2norm(i,j,:)-(nanmean(squeeze(...
                        data2norm(i,j,[frames2norm(1):frames2norm(2),...
                        frames2norm(3):frames2norm(4)])))))...
                        ./abs(nanmean(squeeze(data2norm(i,j,[frames2norm(1):...
                        frames2norm(2),frames2norm(3):frames2norm(4)])))); 
                end %end of check if have multiple norm periods
            end %end dim 2 loop
        end %end dim 1 loop
    end
    %loop through cell, and normalising data within
elseif iscell(data2norm)==1 %inputted data is a cell
    for a=1:size(data2norm,2) %loop cells
        data_ttt=data2norm{1,a}; %output data from within cell (so now mat)
        if size(isnan(data_ttt),1)>1 %check not a nan
            %find out if the data has 2 or 3 dimensions
            t=ndims(data_ttt);
            %again check time dimension is last for 2D data
            if t==2
                if size(data_ttt,1)>size(data_ttt,2)
                    %if not, transpose
                    data_ttt=data_ttt';
                end
            end
            %for 2D data, just need to loop through 1st dim, and norm
            if t==2
                for i=1:size(data2norm,1) %loop dim 1
                    if size(frames2norm,2)==2 %sent in just one norm period
                        %put data back into cell
                        dataNorm{1,a}(i,:)= (data_ttt(i,:)-(nanmean(squeeze(...
                            data_ttt(i,frames2norm(1):frames2norm(2))))))...
                            ./abs(nanmean(squeeze(data_ttt(i,frames2norm(1):...
                            frames2norm(2)))));
                    else
                        dataNorm{1,a}(i,:)= (data_ttt(i,:)-(nanmean(squeeze(...
                            data_ttt(i,[frames2norm(1):frames2norm(2),...
                            frames2norm(3):frames2norm(4)])))))...
                            ./abs(nanmean(squeeze(data_ttt(i,[frames2norm(1):...
                            frames2norm(2),frames2norm(3):frames2norm(4)]))));
                    end
                end %end of loop dim 1
                %for 3D data, need to loop through 1st 2 dims and norm
            elseif t==3
                %normalise data by dividing by the mean of the norm period
                for i=1:size(data_ttt,1) %loop dim 1
                    for j=1:size(data_ttt,2) %loop dim 2
                        if size(frames2norm,2)==2 %sent in just one norm period
                            %put data back into cell
                            dataNorm{1,a}(i,j,:)= (data_ttt(i,j,:)-(nanmean(...
                                squeeze(data_ttt(i,j,frames2norm(1):...
                                frames2norm(2))))))./abs(nanmean(squeeze(...
                                data_ttt(i,j,frames2norm(1):frames2norm(2)))));
                        else
                            %put data back into cell
                            dataNorm{1,a}(i,j,:)= (data_ttt(i,j,:)-(nanmean(...
                                squeeze(data_ttt(i,j,[frames2norm(1):...
                                frames2norm(2),frames2norm(3):...
                                frames2norm(4)])))))./abs(nanmean(squeeze(...
                                data_ttt(i,j,[frames2norm(1):frames2norm(2),...
                                frames2norm(3):frames2norm(4)]))));
                        end
                    end %end of loop dim 1
                end %end of loop dim 2
            end %end of checking how many dimensions in mat within cell
        else %NaN check
            %if data is a NaN, just output a NaN
            dataNorm{1,a}=NaN;
        end %end of checking if cell just has a NaN inside
    end %end of looping cells
end %end of cell check

end %end of function

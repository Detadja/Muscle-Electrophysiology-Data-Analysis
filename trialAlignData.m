function [trial_data] = trialAlignData(data, align_times, time_before, time_after, Fs)

%%%%% All lines where you have to fill in information is tagged with a comment including "FILLIN". Use this flag to find everything you need to modify.
% The function description below described the high-level goals of the function and formats of the inputs and outputs. Read this carefully.


%trialAlignData   (A. Orsborn, created: 4-22-10, updated 1/17/20)
%[trial_data] = trialAlignData(data, align_times, time_before, time_after, Fs)
%aligns a data time-series (data) to the align_times. Returns
%time-series data segments that are [time_before+time_after] long for each trial that
%begin time_before prior to the trial-event
%
%input:
%       data         - data to trial sort. [time x n] matrix, where n is # variables.
%       align_times  - time of trial events to align to [#trials x 1]. Time units must be
%                      consistent with sampling rate (e.g. if event-times in seconds, Fs must be in Hz)
%       time_before  - time before event-code occurence to include in data (in units consistent with Fs)
%       time_after   - total length of trial-aligned data segments (in units consistent with Fs)
%       Fs           - sampling rate
%output:
%       trial_data - Trial-aligned data [trials x time x n]


num_trials  = length(align_times); %# of trials
num_channels = size(data,2);
time_stamps = (time_before + time_after) * Fs;

%initialize trial_data to a [#trials #time-stamps #channels] matrix
%hint: #time-stamps will be the total length of time of a trial
%=(time_before + time_after) converted from units of time to units of samples
trial_data = zeros(num_trials, time_stamps, num_channels); %FILLIN


%loop through trials
for iT=1:num_trials

    %select the indices for this trial
    %note that the 'floor' operation is used to assure we keep the indices
    %as integer values.
    idx_start = floor( (align_times(iT) - time_before)* Fs) ;
    idx_end   = idx_start + (time_before + time_after)*Fs - 1;


    %confirm within bounds of data (only happens for incorrectly formatted data OR
    %when trials occur very close to the start/end of a file.
    if idx_start > 0 & idx_end > 0 & idx_end <= size(data,1) %#ok<AND2>

        %Save the data between your start and end indices into trial_data
        %for the current trial
        trial_data(iT,:,:) = data(idx_start:idx_end, :); %FILLIN

    else
        disp(strcat('Error in time indexing. Skipping trial ', num2str(iT)))
    end
end

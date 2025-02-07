%%%%% This script provides a guided outline for analyzing your experimental data collected
%%%%% for Experiment 3 (characterizing evoked muscle twitches in response to stimulation at differing frequencies)
%%%%% Written by A.L. Orsborn, v210122
%%%%%
%%%%%
%%%%% All lines where you have to fill in information is tagged with a comment including "FILLIN". Use this flag to find everything you need to modify.


%%%%IMPORTANT:
% This will use all the same data processing techniques and calculations we developed in experiment 1.
%Do this analysis AFTER having done the analysis for that experiment.
%you should be able to directly lift most of your code from lab 1 with slight modifications for this analysis.


%% load the data for all files, compute metrics, and make plots


%start with a fresh workspace
clear all

%again, we want to point to where our data files are. Except now, we want to specify a list of all files we want to analyze.
dataDir = ''; %FILLIN with the path to where your data is stored

file_prefix = 'Lab3.2.'; %FILLIN with the text string that is common among all your data files
file_type = '.wav';   %FILLIN with the file extension for your data type

file_date   = ''; %FILLIN with the date string used in your file
file_idnums  = [2,3,4,5,6] ; %FILLIN with the LIST of numeric values of the file ids.
%excluded 1st trial


%we also want to define the meta-data associated with each file we listed so we can analyze trends.
stim_frequency = [3,4,6,8,10]; %FILLIN %excluded 1st trial


num_files = length(file_idnums);


% loop through your files to calculate your evoked twitch profiles and metrics for each.
% here, we will only look at maximum amplitude of the twitch and the absolute force as a function of frequency
mean_twitch_amp = nan(num_files,1);
ste_twitch_amp = nan(num_files,1);
mean_force_amp = nan(num_files,1);
ste_force_amp = nan(num_files,1);
for iF=1:num_files

    %1. load data
    %FILLIN
    [audio_data, sampling_rate] = audioread(fullfile(dataDir, [file_prefix, num2str(file_idnums(iF)), file_type]));  % FILLIN
    force_data = audio_data;  % Replace with actual variable name in your file
    

    %correct the force-sensor data offset
    %FILLIN
    force_offset = mean(force_data(1:100));  % Assuming the first 100 points are baseline
    force_data_corrected = force_data - force_offset;


    %2. detect stim events
    %FILLIN
    stim_threshold = 0.5;  % FILLIN with appropriate threshold based on your data
    stim_events = find(force_data_corrected > stim_threshold);

    %3. trial-align data to stimulation times. Remember to correct for the baseline force in each trial & create a trial_time vector.
    %FILLIN
    trial_time = (1:length(force_data_corrected)) / sampling_rate;  % Adjust for actual sampling rate
    stim_time = trial_time(stim_events);  % Find times of stimulation
    trial_force_corrected = force_data_corrected(stim_events);


    %4. Now compute twitch amplitude and force amplitude for each trial.
    %since we're only calculating the max for each trial, we no longer need a for-loop because the max function can operate across a matrix.
    %look at the help for 'max' to see how to do this calculation without a for-loop
    num_trials = length(stim_events);

    %compute the contraction amplitude (maximum change in force evoked by stim)
    %Here, we want baseline variation removed to see the delta caused by stim.
    %ta = max(); %%FILLIN - twitch amplitude
    ta = max(trial_force_corrected - force_offset); % Twitch amplitude

    %compute the maximum absolute force for the trial
    %this requires using our "uncorrected" force profiles (where we haven't removed the baseline force)
    %fa = max(); %%FILLIN - force amplitude
    fa = max(abs(trial_force_corrected)); % Force amplitude


    %5. now compute trial-averaged metrics to summarize results
    %ste = standard error = standard deviation /sqrt(# trials)
    %hint: look up help for std
    ste_twitch_amp(iF) = std(ta) / sqrt(num_trials);  % Standard error for twitch amplitude
    ste_force_amp(iF) = std(fa) / sqrt(num_trials);  % Standard error for force amplitude

    mean_twitch_amp(iF) = mean(ta);  % Mean twitch amplitude
    mean_force_amp(iF) = mean(fa);  % Mean force amplitude


end %end loop through files.


% generate figures to visualize the relationship between stimulation amplitude and
% each twitch paramter
%plot the across-trial mean with error bars showing the standard error
%hint: look up the help for 'errorbar'
%INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS
figure
subplot(1,2,1)
%errorbar() %FILLIN
errorbar(stim_frequency, mean_twitch_amp, ste_twitch_amp, 'o')
xlabel('Stimulation frequency (Hz)')
ylabel('Twitch amplitude (V)')
title('Twitch Amplitude vs. Stimulation Frequency')
grid on


subplot(1,2,2)
%errorbar() %FILLIN
errorbar(stim_frequency, mean_force_amp, ste_force_amp, 'o')
xlabel('Stimulation frequency (Hz)')
ylabel('Force amplitude (V)')
title('Force Amplitude vs. Stimulation Frequency')
grid on

%%%%% This script provides a guided outline for analyzing your experimental data collected
%%%%% for Experiment 1 (characterizing evoked muscle twitches in response to stimulation at differing amplitudes)
%%%%% Written by A.L. Orsborn, v210122
%%%%%
%%%%%
%%%%% All lines where you have to fill in information is tagged with a comment including "FILLIN". Use this flag to find everything you need to modify.



% we will first load a data file and test our pre-processing and calculations on one file. Once that is complete, we can extend our analysis to all data to examine trends.
%% testing pre-processing and calculations

% Define some basic things to make it easy to find your data files.
% We will want to take advantage of systematic naming structure in our data files.
% Your files should have names like [prefix][date][id #].
% Note that the SpikeRecorder program automatically saves files with date and time in the name.
% We recommend re-naming your files to convert time into a simpler id# e.g. 1, 2, 3...

dataDir = 'C:\Users\trnan\OneDrive\Desktop\uw\EE\EE 466\Lab3\Experimental Data\Experiment 1\'; %%FILLIN with the path to where your data is stored

file_prefix = 'Lab3.1'; %%FILLIN with the text string that is common among all your data files
file_type = '.wav';   %%FILLIN with the file extension for your data type

file_date   = ''; %%FILLIN with the date string used in your file
file_idnum  = '3'; %%FILLIN with the numeric value of the file id

full_file_name = [dataDir file_prefix file_date '.' num2str(file_idnum) file_type];

% load your data file and the sampling rate into matlab into the variables 'data' and FS, respectively.
[data, FS] = audioread(full_file_name);

% check the basic properties of your loaded variables
whos data FS

% plot 20 seconds of your data to look at the data scales
figure
plot( data([1: FS * 20],: ) ) %FILLIN the portion in the brackets to only plot 20 seconds
xlabel('Time, in samples')
ylabel('Voltage (mV)')

%amplifier settings in the spikerbox analog inputs add a constant offset voltage, which makes our forces appear negative.
%we will use the minimum observed voltage on the analog channel (channel 2) to shift it to be > 0
force_voltage_offset = min(data(:,2));
data(:,2) = data(:,2) - force_voltage_offset; %FILLIN -- use force_voltage_offset to adjust data(:,2) so that it is > 0

%now re-plot to make sure the force trace is now > 0
figure
plot( data )
xlabel('Time, in samples')
ylabel('Voltage (mV)')



%now we need to use our channel 1 inputs to detect stimulation events.
%

threshold = 0.8; %FILLIN--you may need to modify the threshold for your data--look at previous plots to determine the best threshold to use
[stim_times, stim_waveforms] = detectSpikes(data(:,1), threshold, FS);

%2. The stimulation waveform is complex and will cross our threshold several times for a single stimulation event.
%To remove these 'extra' detected stimulation times, we can use our knowledge that stimulation pulse timing.
%We know that stimulation pulses should not be _very_ close together.
MIN_STIM_DELTA = 0.2; % FILLIN - make sure the MIN_STIM_DELTA is less than the actual Stim freq (e.g. 0.3 < 1 / 2Hz)
dt_stim = diff(stim_times); %use 'diff' to find the time between detected stimulation pulses
too_short = find(dt_stim<MIN_STIM_DELTA); %find pulses that are too close together.
use_pulse = setdiff(1:length(stim_times), too_short); %pulses to keep are all the other ones. setdiff(A,B) returns all the values in the vector A that are NOT in vector B.
stim_times_cleaned = stim_times(use_pulse);


%now that we have the stimulation times, we can use them to trial-align the force pulses
%to the time of stimulation. 
%This operation (logic and steps) are the same as when we trial-aligned
%spike times, but we now have a continuous variable (a force)
%rather than a discrete variable (spike times).
%We've provided a shell function trialAlignData.m to perform these
%computations. Open it and finish filling it in.
align_times = stim_times_cleaned(:,1); %FILLIN
time_before = 0.05;
time_after  = 0.3; %the twitches are very short--we don't need to load much time.
[trial_force] = trialAlignData(data(:,2), align_times, time_before, time_after, FS); %FILLIN. reminder: the force data is the second column in your data set. We only want to align the force data.

%create a time-axis for our trial-aligned data (since trial_data itself is
%in samples)
trial_time = linspace(-time_before, time_after, size(trial_force,2));


%plot the trial-aligned data to check our code worked properly.
figure
plot(trial_time, trial_force)
xlabel('Time (s)')
ylabel('Force sensor voltage, (V)')
hold on
title('Trial-aligned twitch force (raw)')


%In your plot above, it's clear there's a shift in the average force across trials.
%This does not come from the stimulation, but from other sources of variability.
%for instance, how much the subject is pushing on the sensor may vary across time, shifting the 'baseline' force
%The twitch force profile is observed on top of that drifting baseline.
%Let's clean up our force pulse data a bit to remove this extra variability in baseline force.
%We will do this by subtracting the force at the time of stimulation (i.e. the 'baseline' force) from each trial individually
%This zeros out any differences in baseline across trials
baseline_idx = find(trial_time>=0, 1, 'first');
trial_baseline_force = trial_force(:, baseline_idx); %this is [#trials x 1] vector
trial_force_corrected = trial_force - repmat(trial_baseline_force, [1,1750]); %FILLIN to subtract the baseline force for each trial from all time-points in the trial.
%hint: look at the help for 'repmat'. Use it to duplicate trial_baseline_force in the right way so that it's [#trials #time_points].

%now plot the trial-aligned data to check our code worked properly.
%INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS
figure
plot(trial_time, trial_force_corrected)
xlabel('Time (s)')
ylabel('Force sensor voltage, (V)')
hold on
plot([0 0], [0 max(trial_force_corrected(:))], 'k--') %plots a dashed line at t=0 to visualize alignment
title('Trial-aligned twitch force, baseline corrected')



%3. Now let's write code to compute the relevant metrics of twitch force profiles
%we will compute these metrics on each individual trial and then average across trials.

num_trials = size(trial_force_corrected,1);

%initialize vectors for metrics
latent_period = nan(num_trials,1);
contraction_time = nan(num_trials,1);
relaxation_time = nan(num_trials,1);
contraction_amp = nan(num_trials,1);

%loop through trials
for iT=1:num_trials

  %compute the latency. We will define this by computing a pre-stimulation 'baseline'
  %and finding the first time when the force > 5% baseline
  baseline = mean(abs(trial_force_corrected(iT, trial_time<0)),2); %make sure you understand what this is doing.

  %Find the first point when force exceeds threshold. Only looking at data AFTER stimulation
  idx_rise = find(trial_force_corrected(iT, trial_time>0) > 0.05 * baseline, 1, 'first'); %FILLIN to find when force exceeds our specified threshold.
  idx_rise = idx_rise + find(trial_time>0, 1, 'first')-1; %adjust idx_rise to match indexing of full data matrix (since we searched a subset)

  %to make code robust to noisy data, check that you find an onset time. Is not guaranteed.
  if ~isempty(idx_rise)
      latent_period(iT) = trial_time(idx_rise);
  else
    warning('Cannot find latent period')
    continue
    %if can't define latent period, the value of latent_period(iT) will remain a nan (missing data)
  end

  %compute the contraction amplitude (maximum force evoked) and contraction time (time it takes to hit max)
  [mx, idx_mx] = max(trial_force_corrected(iT,:)); %FILLIN

  contraction_amp(iT) = mx;
  contraction_time(iT) = (idx_mx - idx_rise) * 1 / FS; %FILLIN. Be sure this is in units of time, not samples.

  %compute the relaxation_time as the point when force falls below 60% of maximum
  %note, have to make sure this point is AFTER maximum contraction time
  %Use the same logic/technique as what we did above for idx_rise.
  idx_rt = find(trial_force_corrected(iT, [1:1750] > idx_mx) < contraction_amp(iT)*.6,1, 'first'); %FILLIN the brackets to search the right piece of the data
  idx_rt = idx_rt + idx_mx; %FILLIN to account for searching a subset of the data.

  %to make code robust to noisy data, check that you find a relaxation time. Is not guaranteed.
  if ~isempty(idx_rt)
      relaxation_time(iT) = trial_time(idx_rt);
  else
    warning('Cannot find relaxation time')
    continue
    %if isempty(idx_rt), relaxation_time for that trial will remain a nan, i.e. missing data
  end

end %end loop of trials.

%display trial-averaged results
disp(['Mean contraction amplitude:', num2str(nanmean(contraction_amp)), 'volts'])
disp(['Mean contraction time:', num2str(nanmean(contraction_time)*1000), 'milliseconds']) %note: convert from seconds to ms
disp(['Mean relaxation time:', num2str(nanmean(relaxation_time)*1000), 'milliseconds'])
disp(['Mean latent period:', num2str(nanmean(latent_period)*1000), 'milliseconds'])


%also plot trial-averaged twitch profile (average across trials of trial_force_corrected)
%INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS
figure
plot(trial_time, mean(trial_force_corrected, 1))%FILLIN
hold on
plot([0 0], [0 max(trial_force_corrected(:))], 'k--') %plots a dashed line at t=0 to visualize alignment.
%Fill in the brackets with the data vector you're plotting.
xlabel('Time (seconds)')
ylabel('Twitch force (V)')

%% Now we will load the data for all files, compute metrics, and make plots


%start with a fresh workspace
clear all

%again, we want to point to where our data files are. Except now, we want to specify a list of all files we want to analyze.
dataDir = 'C:\Users\trnan\OneDrive\Desktop\uw\EE\EE 466\Lab3\Experimental Data\Experiment 1\'; %FILLIN with the path to where your data is stored

file_prefix = 'Lab3.1'; %FILLIN with the text string that is common among all your data files
file_type = '.wav';   %FILLIN with the file extension for your data type

file_date   = ''; %FILLIN with the date string used in your file
file_idnums  = [1 2 3 4] ; %FILLIN with the LIST of numeric values of the file ids.


%we also want to define the meta-data associated with each file we listed so we can analyze trends.
stim_amplitude = [4 4.25 4.5 4.75]; %FILLIN


num_files = length(file_idnums);


% loop through your files to calculate your evoked twitch profiles and metrics for each.

%initialize variables for mean and standard-error metrics across files
%ste = standard error = standard deviation /sqrt(# trials)
mean_contraction_amp = nan(num_files,1);
ste_contraction_amp = nan(num_files,1);
mean_contraction_time = nan(num_files,1);
ste_contraction_time  = nan(num_files,1);
mean_relaxation_time = nan(num_files,1);
ste_relaxation_time = nan(num_files,1);
mean_latent_period = nan(num_files,1);
ste_latent_period = nan(num_files,1);
std_error = nan(num_files,1);

for iF=1:num_files

    %For these steps, you can directly repurpose code you wrote above when testing a single file.
    %Only modify to change the file you're loading.
    full_file_name = [dataDir file_prefix file_date '.' num2str(file_idnums(iF)) file_type];

    %1. load data
    %FILLIN
    [data, FS] = audioread(full_file_name);

    %correct the force-sensor data offset
    %FILLIN
    force_voltage_offset = min(data(:,2));
    data(:,2) = data(:,2) - force_voltage_offset; %FILLIN -- use force_voltage_offset to adjust data(:,2) so that it is > 0

    %2. detect stim events
    %FILLIN
    threshold = 0.8; %FILLIN--you may need to modify the threshold for your data--look at previous plots to determine the best threshold to use
    [stim_times, stim_waveforms] = detectSpikes(data(:,1), threshold, FS);
    MIN_STIM_DELTA = 0.3; % FILLIN - make sure the MIN_STIM_DELTA is less than the actual Stim freq (e.g. 0.3 < 1 / 2Hz)
    dt_stim = diff(stim_times); %use 'diff' to find the time between detected stimulation pulses
    too_short = find(dt_stim<MIN_STIM_DELTA); %find pulses that are too close together.
    use_pulse = setdiff(1:length(stim_times), too_short); %pulses to keep are all the other ones. setdiff(A,B) returns all the values in the vector A that are NOT in vector B.
    stim_times_cleaned = stim_times(use_pulse);

    %3. trial-align data to stimulation times. Remember to correct for the baseline force in each trial & create a trial_time vector.
    %FILLIN
    align_times = stim_times_cleaned(:,1); %FILLIN
    time_before = 0.05;
    time_after  = 0.3; %the twitches are very short--we don't need to load much time.
    [trial_force] = trialAlignData(data(:,2), align_times, time_before, time_after, FS); %FILLIN. reminder: the force data is the second column in your data set. We only want to align the force data.
    
    %create a time-axis for our trial-aligned data (since trial_data itself is
    %in samples)
    trial_time = linspace(-time_before, time_after, size(trial_force,2));
    
    %In your plot above, it's clear there's a shift in the average force across trials.
    %This does not come from the stimulation, but from other sources of variability.
    %for instance, how much the subject is pushing on the sensor may vary across time, shifting the 'baseline' force
    %The twitch force profile is observed on top of that drifting baseline.
    %Let's clean up our force pulse data a bit to remove this extra variability in baseline force.
    %We will do this by subtracting the force at the time of stimulation (i.e. the 'baseline' force) from each trial individually
    %This zeros out any differences in baseline across trials
    baseline_idx = find(trial_time>=0, 1, 'first');
    trial_baseline_force = trial_force(:, baseline_idx); %this is [#trials x 1] vector
    trial_force_corrected = trial_force - repmat(trial_baseline_force, [1,1750]); %FILLIN to subtract the baseline force for each trial from all time-points in the trial.
    
    %4. Now compute twitch metrics for each trial.
    num_trials = size(trial_force_corrected,1);

    %initialize vectors for metrics
    lp = nan(num_trials,1);
    ct = nan(num_trials,1);
    rt = nan(num_trials,1);
    ca = nan(num_trials,1);

    %loop through trials
    for iT=1:num_trials

      %compute the latency. We will define this by computing a pre-stimulation 'baseline'
      %and finding the first time when the force > 5% baseline
      baseline = mean(abs(trial_force_corrected(iT, trial_time<0)),2); %make sure you understand what this is doing.
    
      %Find the first point when force exceeds threshold. Only looking at data AFTER stimulation
      idx_rise = find(trial_force_corrected(iT, trial_time>0) > 0.05 * baseline, 1, 'first'); %FILLIN to find when force exceeds our specified threshold.
      idx_rise = idx_rise + find(trial_time>0, 1, 'first')-1; %adjust idx_rise to match indexing of full data matrix (since we searched a subset)
    
      %to make code robust to noisy data, check that you find an onset time. Is not guaranteed.
      if ~isempty(idx_rise)
          lp(iT) = trial_time(idx_rise); %FILLIN
      else
        warning('Cannot find latent period')
        continue
        %if can't define latent period, the value of latent_period(iT) will remain a nan (missing data)
      end

      %compute the contraction amplitude (maximum force evoked) and contraction time (time it takes to hit max)
      %compute the contraction amplitude (maximum force evoked) and contraction time (time it takes to hit max)
      [mx, idx_mx] = max(trial_force_corrected(iT,:)); %FILLIN
    
      ca(iT) = mx;  %FILLIN
      ct(iT) = (idx_mx - idx_rise) * 1 / FS; %FILLIN. Be sure this is in units of time, not samples.

      %compute the relaxation_time as the point when force falls below 60% of maximum
      %note, have to make sure this point is AFTER maximum contraction time
      idx_rt = find(trial_force_corrected(iT, [1:1750] > idx_mx) < ca(iT)*.6,1, 'first'); %FILLIN the brackets to search the right piece of the data
      idx_rt = idx_rt + idx_mx; %FILLIN to account for searching a subset of the data.
      %to make code robust to noisy data, check that you find a relaxation time. Is not guaranteed.
      if ~isempty(idx_rt)
          rt(iT) = trial_time(idx_rt); %FILLIN
      else
        warning('Cannot find relaxation time')
        continue
        %if isempty(idx_rt), relaxation_time for that trial will remain a nan, i.e. missing data
      end

    end %end loop of trials.


    %5. now compute trial-averaged metrics to summarize results
    %ste = standard error = standard deviation /sqrt(# trials)
    %hint: look up help for nanstd
    mean_contraction_amp(iF) = nanmean(ca);
    ste_contraction_amp(iF) = nanstd(ca) / sqrt(sum(~isnan(ca))); %FILLIN

    mean_contraction_time(iF) = nanmean(ct);
    ste_contraction_time(iF)  = nanstd(ct) / sqrt(sum(~isnan(ct))); %FILLIN

    mean_relaxation_time(iF) = nanmean(rt);
    ste_relaxation_time(iF) = nanstd(rt) / sqrt(sum(~isnan(rt))); %FILLIN

    mean_latent_period(iF) = nanmean(lp);
    ste_latent_period(iF) = nanstd(lp) / sqrt(sum(~isnan(lp))); %FILLIN


    %also save the mean twitch trace to plot overlayed
    mean_twitch_profile(iF,:) = mean(trial_force_corrected, 1); %FILLIN
    std_error(iF) = nanstd(mean_twitch_profile(iF, :)) ./ sqrt(size(mean_twitch_profile, 2)); % Standard error

end %end loop through files.


% generate figures to visualize the relationship between stimulation amplitude and
% each twitch paramter
%plot the across-trial mean with error bars showing the standard error
%hint: look up the help for 'errorbar'
%INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS
figure
subplot(1,4,1)
errorbar(mean_contraction_amp, ste_contraction_amp) %FILLIN
% plot(trial_time,mean_twitch_profile(1, :))%FILLIN
xlabel('Stimulation amplitude (a.u.)')
ylabel('Contraction amplitude (V)')

subplot(1,4,2)
errorbar(mean_contraction_time, ste_contraction_time) %FILLIN
% plot(trial_time, mean_twitch_profile(2, :))%FILLIN
xlabel('Stimulation amplitude (a.u.)')
ylabel('contraction time (ms)')

subplot(1,4,3)
errorbar(mean_relaxation_time, ste_relaxation_time) %FILLIN
% plot(trial_time, mean_twitch_profile(3, :))%FILLIN
xlabel('Stimulation amplitude (a.u.)')
ylabel('Relaxation time (ms)')

subplot(1,4,4)
errorbar(mean_latent_period, ste_latent_period) %FILLIN
% plot(trial_time, mean_twitch_profile(4, :))%FILLIN
xlabel('Stimulation amplitude (a.u.)')
ylabel('Latent period (ms)')


%finally, plot all the twitches overlayed
%INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS
figure
plot(trial_time, mean_twitch_profile)
xlabel('Time (s)')
ylabel('Twitch Force (V)')
hold on
plot([0 0], [min(mean_twitch_profile(:)) max(mean_twitch_profile(:))], 'k--')
legend({'4.00 V' '4.25 V' '4.50 V' '4.75 V' 'onset'}) %FILLIN

%%%%% This script provides a guided outline for analyzing your experimental data collected
%%%%% for Experiment 1 (quantifying EMG's relationship to force)
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

dataDir = 'C:\Users\denis\OneDrive\Uni Stuff\UW\2 Winter 2024\BIOEN 566 - Neural Engineering Lab\Lab Manuals and Code Templates\Lab 3\Experimental Data\Lab 3 Exp 3 Data - Luke\'; %FILLIN with the path to where your data is stored
% dataDir = 'C:\Users\denis\OneDrive\Uni Stuff\UW\2 Winter 2024\BIOEN 566 - Neural Engineering Lab\Lab Manuals and Code Templates\Lab 3\Experimental Data\Experiment 3\';

file_prefix = 'BYB_Recording_2021-01-27_'; %FILLIN with the text string that is common among all your data files
% file_prefix = 'Lab3.3.';
file_type = '.wav';   %FILLIN with the file extension for your data type

% file_date   = ''; %FILLIN with the date string used in your file
file_idnum  = 11; %FILLIN with the numeric value of the file id

full_file_name = [dataDir file_prefix num2str(file_idnum) file_type];

% load your data file and the sampling rate into matlab into the variables 'data' and FS, respectively.
% hint: look at the Matlab function 'audioread'
[data, FS] = audioread(full_file_name);

% check the basic properties of your loaded variables
whos data FS


% plot 20 seconds of your data to assure the loaded file looks correct
figure
plot( data ) %FILLIN the portion in the brackets to only plot 20 second
xlabel('Time, in samples')
ylabel('Voltage (mV)')


%now let's trial-align our EMG to the time of force onset

%1. load an 'events' file.
%Use the readEventsFile.m function you used in lab 2 to load up the data
events_file_type = '.txt';
events_file_name = [dataDir file_prefix num2str(file_idnum) events_file_type];
[EVENTS, EVENT_TIMES] = readEventsFile(events_file_name);


%2. Align your data to the trial-events.
%you already wrote the trialAlignData function for experiment 1 analysis.
%We can use that function here too, because EMG is also a continuous
%signal like the force. 
align_times = EVENT_TIMES;
time_before = 1;
time_after  = 2; %set this based on how long you held your contractions
[trial_data] = trialAlignData(data, align_times, time_before, time_after, FS);

%create a time-axis for our trial-aligned data (since trial_data itself is
%in samples)
trial_time_emg = linspace(-time_before, time_after, size(trial_data,2));


%now plot the trial-aligned data to check our code worked properly.
%INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS
figure
plot(trial_time_emg, trial_data');
hold on
plot([0 0], [0 max(trial_data(:))], 'k--') %plots a dashed line at t=0 to visualize alignment
xlabel('Time (s)')
ylabel('EMG (Volts)')
title('Trial-aligned EMG')



%3. Now let's write code to compute the relevant metrics of EMG that we
%will test for our analyses

%3a - mean EMG in a time-window
avg_start = 0; %time to start (in seconds)
avg_end   = 5; %time to end averaging (in seconds) - again, set this based on how long you held your contractions
avg_idx   = trial_time_emg >= avg_start & trial_time_emg <= avg_end; %FILLIN - create a logical vector to indicate the indices of trial_data where trial_time_emg is between your avg_start and avg_end times.

%Hint: mean can operate on a matrix to compute the mean over 1 dimension.
%You want to average over time in your [trials x time] matrix. Use avg_idx created above to
%only take the average over a defined subset of time
emg_amp = mean(trial_data(:, avg_idx), 2 ); %FILLIN.


%plot a histogram of the mean emg amplitude across your trials to visualize results
figure
hist(emg_amp,10)
xlabel('Mean EMG Amplitude') %FILLIN
ylabel('Frequency of Amplitudes') %FILLIN
title('Distribution of EMG ampiltude across trials')



%3b - compute the RMS value (root-mean-square) as a measure of EMG signal variance
%in a given time-window.
%hint: type help for the matlab function 'rms'. It can be used the same way you used 'mean' above.
rms_start = 0; %time to start (in seconds)
rms_end   = 5; %time to end averaging (in seconds) - again, set this based on how long you held your contractions
rms_idx   = trial_time_emg >= rms_start & trial_time_emg <= rms_end; %FILLIN - create a logical vector to indicate the indices of trial_data where trial_time_emg is between your rms_start and rms_end times.

emg_rms = rms(trial_data(:, rms_idx), 2); %FILLIN

%plot a histogram of the RMS
figure
hist(emg_rms,10)
xlabel('EMG RMS Amplitudes') %FILLIN
ylabel('Frequency of Amplitudes') %FILLIN
title('Distribution of EMG RMS across trials')


%3c - compute the spectral power of the EMG signal
%we have provided a COMPLETED function to compute the periodogram (spectral
%power estimate) for each trial called 'compute_periodogram_fft.m'
%look at the help for this function and use it to analyze your data

%pick time-window for analysis (once force has started)
spec_start = 0; %time to start (in seconds)
spec_end   = 5; %time to end averaging (in seconds) - again, set this based on how long you held your contractions
spec_idx   =  trial_time_emg >= spec_start & trial_time_emg <= spec_end; %FILLIN - create a logical vector to indicate the indices of trial_data where trial_time_emg is between your spec_start and spec_end times.

%note: the muscleBox includes hardware filtering to reduce low-frequency activity.
%We therefore recommend looking at frequencies >70 hz and less than 500 Hz.
fpass = [70 500];
[emg_periodogram, spec_frequency] = compute_periodogram_fft(trial_data(:, spec_idx), FS, fpass); %FILLIN. Remember to only analyze the defined time-window from spec_idx.


%plot the average spectral power across trials
%plot power vs. frequency to visualize
figure
plot(spec_frequency, mean(emg_periodogram,1))
xlabel('Frequency (Hz)')
ylabel('EMG power (mV^2)')


%% Now we will load the data for all files, compute metrics, and make plots


%start with a fresh workspace
clear all

%again, we want to point to where our data files are. Except now, we want to specify a list of all files we want to analyze.
%dataDir = 'C:\Users\Orsborn Lab\Desktop\Amy\lab3\data\200118\'; %FILLIN with the path to where your data is stored
dataDir = 'C:\Users\denis\OneDrive\Uni Stuff\UW\2 Winter 2024\BIOEN 566 - Neural Engineering Lab\Lab Manuals and Code Templates\Lab 3\Experimental Data\Lab 3 Exp 3 Data - Luke\'; %FILLIN with the path to where your data is stored

file_prefix = 'BYB_Recording_2021-01-27_'; %FILLIN with the text string that is common among all your data files
file_type = '.wav';   %FILLIN with the file extension for your data type
events_file_type = '.txt'; %file-flag + extension for event data

% file_date   = ''; %FILLIN with the date string used in your file
file_idnums  = [11 12 13 14] ; %FILLIN with the LIST of numeric values of the file ids.


%we also want to define the meta-data associated with each file we listed so we can analyze trends.
%for reasons related to later analysis, be sure this is a ROW vector (i.e.
%size = [#forces x 1]
grip_force = [10 20 30 40]; %FILLIN

%for spectral analysis, we'll analyze 3 frequency bands:
fbands(1,:) = [100 200];
fbands(2,:) = [200 300];
fbands(3,:) = [300 400];
num_bands = size(fbands,1);

%initialize matrices for putting together data across files. ---NOTE: DON'T
%NEED TO DO THIS BECAUSE CONCAENATING
num_files = length(file_idnums);
% emg_amps = []; %FILLIN
% emg_rmss = []; %FILLIN
% emg_powers = cell(num_files,1);
% trial_file_id = []; %%FILLIN - this will be an ID to keep track of which file each data point belongs to.

%initialize variables for mean and standard-error metrics across files
%ste = standard error = standard deviation /sqrt(# trials)
emg_mean_amp = nan(num_files,1);
emg_mean_rms = nan(num_files,1);
emg_mean_power = nan(num_files, num_bands);
emg_ste_amp = nan(num_files,1);
emg_ste_rms = nan(num_files,1);
emg_ste_power = nan(num_files, num_bands);

% loop through your files to calculate your EMG metrics for each grip-force level
for iF=1:num_files

    %load emg data
    %FILLIN
    full_file_name = [dataDir file_prefix num2str(file_idnums(iF)) file_type];
    [data, FS] = audioread(full_file_name);

    %load events
    %FILLIN
    events_file_name = [dataDir file_prefix num2str(file_idnums(iF)) events_file_type];
    [EVENTS, EVENT_TIMES] = readEventsFile(events_file_name);

    %trial-align
    %FILLIN
    align_times = EVENT_TIMES;
    time_before = 1;
    time_after  = 2;
    [trial_data] = trialAlignData(data, align_times, time_before, time_after, FS);

    %create a time-axis for our trial-aligned data
    %FILLIN
    trial_time_emg = linspace(-time_before, time_after, size(trial_data,2));

    %compute mean (across time for each trial--like we did above.)
    %you can re-purpose the code from above here.
    avg_start = 0; 
    avg_end   = 5;
    avg_idx   = trial_time_emg >= avg_start & trial_time_emg <= avg_end;
    amp_tmp = mean( trial_data(:, avg_idx), 2 ); %%FILLIN - mean in time-window for each trial. Should be [# trials x 1] vector

    %concatonate current file's data with all other data across files.
    % emg_mean_amp = cat(1, emg_mean_amp, amp_tmp);

    %Also keep track of which file (and therefore grip_force) each trial belongs to in trial_file_id
    %hint: if file 1 has 2 trials and file 2 has 5 trials, trial_file_id = [1 1 2 2 2 2 2]
    % trial_file_id = cat(1, trial_file_id, [1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4]); %FILLIN the brackets with the appropriate vector.

    %compute the mean and standard error across trials (i.e. mean of (amp_tmp) as well.
    %standard error (ste) = standard deviation / sqrt(# data points)
    emg_mean_amp(iF) = mean(amp_tmp); %FILLIN.
    emg_ste_amp(iF) = std( amp_tmp )./sqrt(length(amp_tmp)); %FILLIN


    %Repeat the same computations as above, but to compute the RMS.
    rms_start = 0;
    rms_end   = 5;
    rms_idx   = trial_time_emg >= rms_start & trial_time_emg <= rms_end;
    rms_tmp = rms(trial_data(:, rms_idx), 2); %%FILLIN - RMS of time-window for each trial. Should be [# trials x 1] vector

    %concatonate across trials (like we did above)
    % emg_mean_rms = cat(1, emg_mean_rms, rms_tmp); %FILLIN
    %note: don't need another trial_file_id because it's the same for all metrics

    %comute mean and standard error across trials (like above)
    emg_mean_rms(iF) = mean(rms_tmp); %FILLIN
    emg_ste_rms(iF) = std(rms_tmp)./sqrt(length(rms_tmp)); %FILLIN


    %compute the EMG periodogram for the time-band of interest -- reuse the code above.
    %FILLIN
    spec_start = 0; 
    spec_end   = 5; 
    spec_idx   =  trial_time_emg >= spec_start & trial_time_emg <= spec_end; 
    
    fpass = [70 500];
    [emg_periodogram, spec_frequency] = compute_periodogram_fft(trial_data(:, spec_idx), FS, fpass); 
    


    %loop through each frequency - band of interest to get the average power in each band.
    for iB=1:num_bands

      %compute mean within the frequency band.
      freq_idx = spec_frequency >= fbands(iB, 1) & spec_frequency <= fbands(iB, 2); %%FILLIN - create a logical vector that indicates when spec_frequency is in the range [fbands(iB,1) fbands(iB,2)].
      tmp_pow = mean( emg_periodogram(:,freq_idx),2);

      %initialize emg_powers{iB} to empty first time
      % if iF==1
      %   emg_powers{iB} = []; %FILLIN
      % end
      %concatonate across trials (like we did above)
      % emg_powers{iB} = cat(1, emg_powers{iB}, tmp_pow);

      %compute mean and standard error across trials (like above)
      %now, emg_mean_power has two dimensions--the file # and the frequency band.
      emg_mean_power(iF, iB) = mean(tmp_pow); %FILLIN
      emg_ste_power(iF,iB)  = std(tmp_pow) / sqrt(length(tmp_pow));  %FILLIN

    end %end loop through frequency bands


end %end loop through files.

%%
% generate figures to visualize the relationship between grip-force and
% each EMG feature
%plot the across-trial mean with error bars showing the standard error
%hint: look up the help for 'errorbar'
%INCLUDE THIS FIGURE IN COMPREHENSION QUESTIONS
figure
subplot(1,3,1)
errorbar(grip_force, emg_mean_amp, emg_ste_amp) %FILLIN
xlabel('Grip force (N)')
ylabel('EMG Amplitude (mV)')

subplot(1,3,2)
errorbar(grip_force, emg_mean_rms, emg_ste_rms) %FILLIN
xlabel('Grip force (N)')
ylabel('EMG RMS (mV)')

subplot(1,3,3)
%loop through each frequency band of interest.
for iB=1:num_bands
  errorbar(grip_force, emg_mean_power(:,iB), emg_ste_power(:,iB)) %FILLIN
  hold on
end
xlabel('Grip force (N)')
ylabel('EMG spectral power (mV^2)')
legend('100-200 Hz', '200-300 Hz', '300-400 Hz') %FILLIN


%finally, perform a statistical test to confirm which metric shows the most
%linear relationship to grip-force. Here, you should use the single-trial data
%rather than the across-trial averages.
%hint: use the matlab function 'regress' to fit a linear model between the
%metric and grip-force.
%Look at the help for 'regress'.
%Your Y data will be the EMG metric [#total trials across all files x 1]
%Your X data will be [grip_force for each trial; ones(#total trials,1)]. This is a [#tital trials  x 2] matrix
%we add the constant '1' values to allow our linear model to have a
%non-zero intercept.
%
%Hint for creating the X matrix quickly (i.e. without a for-loop):
%If I have a vector A = [5 4 3 2 1], and a vector B = [1 1 2 2 3 3 4 4]
%I can index vector A by the vector B by calling  "A(B)" which returns [5 5 4 4 3 3 2 2] for the example above.
%You can do the same procedure for grip_force and trial_file_id to make X.
%
%To test how well the linear model fits, use the variance explained by the
%linear model (R^2). This is computed by 'regress'. Look at the help to see
%how to get that value back.

%regression for emg-amp
% X = []; %FILLIN
% Y = ; %FILLIN
% [~, ~, ~, ~, stats] = regress(Y, X);
% r2_amp = stats(1);
% 
% disp(['Corr Coef, mean = ', num2str(r2_amp)])
% 
% 
% %change Y for emg_rms
% Y = ; %FILLIN
% [~, ~, ~, ~, stats] = regress(Y, X);
% r2_rms = stats(1);
% 
% disp(['Corr Coef, rms = ', num2str(r2_rms)])
% 
% 
% %change Y for emg power
% for iB=1:num_bands
%   Y = ; %FILLIN
%   [~, ~, ~, ~, stats] = regress(Y, X);
%   r2_power(iB) = stats(1);
% 
%   disp(['Corr Coef, power band ' num2str(iB)' ' = ', num2str(r2_power(iB))])
% end
%BE SURE TO INCLUDE THESE RESULTS IN YOUR RESULTS WRITE-UP

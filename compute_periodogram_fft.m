function [periodogram, frequency] = compute_periodogram_fft(data, Fs, fpass)

%compute_periodogram_fft   (A. Orsborn, updated:1/17/20)
%
%computes periodogram of continuous data using the matlab fast-fourier transform (fft)
%
%input: data         - data (trials x time)
%       Fs           - sampling rate (default: 1000 Hz)
%       fpass        - frequency band of fft to save [min max] (Hz)
%                      default: [0 samplingrate/2]
%output:
%       periodogram - periodogram (trials x freq)
%       frequency   - frequency vector of spectrogram


%if not defined in call, set default parameters
if ~exist('Fs', 'var')
    Fs = 1000;
end
if ~exist('fpass', 'var')
    fpass = [0 Fs/2];
else %argument checking
    if length(fpass)~=2
        error('fpass must be a (1 x 2) vector')
    end
end



if ~isequal(size(data,3), 1) %input checking
    error('data must be (trials x time)')
end

time_window = size(data,2); %# time samples  
frequency = 0:Fs/time_window:Fs/2;

%take only fpass
f1 = find(frequency>=fpass(1), 1, 'first');
f2 = find(frequency<=fpass(2), 1, 'last');
frequency = frequency(f1:f2);




%compute spectral power using fft
%note this computes spectrum for all trials
spec = fft(data,[], 2);
spec = spec(:,1:time_window/2+1); %keep real portion of frequency axis
spec = abs(spec).^2 * (1/Fs*time_window); %convert to power
spec(2:end -1) = spec(2:end -1)*2;

%save the computed data
%return only fpass
periodogram = spec(:,f1:f2);



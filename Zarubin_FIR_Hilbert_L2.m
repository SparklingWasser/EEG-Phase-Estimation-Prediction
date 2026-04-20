function output = Zarubin_FIR_Hilbert_L2(data, fs, target_freq, cutoff_freq, nSamp_predict)
% data: single-channel data (non-filtered)
% fs: sampling rate
% target_f: target frequency, or IAF
% cutoff_freq: BPF cutoff frequencies
% nSamp_predict: number of samples to be predicted (further after the segment lost by filter delay)
% 
% Written by Seonghun Park (s.park7532@gmail.com) according to the article below:
% Zarubin, G., Gundlach, C., Nikulin, V., Villringer, A., & Bogdan, M. (2020). Transient amplitude modulation of alpha-band oscillations by short-time intermittent closed-loop tACS. Frontiers in human neuroscience, 14, 366.
%
% Every parameter is set up as in the reference

data_len = length(data)/fs;

filt_len_in_s = 0.1;
B = fir1(filt_len_in_s*fs, cutoff_freq/(fs/2));

data_filtered = filter(B, 1, data);   % zero-phase filter from literature
data_filtered_delay_corrected = data_filtered(filt_len_in_s*fs/2+1:end);
t_delay_corrected = 1/fs:1/fs:data_len-filt_len_in_s/2;


phase_train = angle(hilbert(data_filtered_delay_corrected));

phase_shift = linspace(0, 2*pi, 72+1);  % 5 degree increment

sin_waves = zeros(length(phase_shift), length(t_delay_corrected));
L2_dist = zeros(length(phase_shift), 1);
for i = 1:length(phase_shift)
    sin_waves(i, :) = sin(2*pi*target_freq*t_delay_corrected+phase_shift(i));
    L2_dist(i) = sqrt(sum(abs(angle(hilbert(sin_waves(i, :))) - phase_train).^2));
end

[~, min_idx] = min(L2_dist);

t_prediction = data_len+1/fs:1/fs:data_len + nSamp_predict/fs;

output = sin(2*pi*target_freq*t_prediction + phase_shift(min_idx));




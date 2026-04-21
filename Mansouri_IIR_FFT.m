function output = Mansouri_IIR_FFT(data, fs, cutoff_freq, nSamp_predict)
% data: single-channel data (non-filtered)
% fs: sampling rate
% cutoff_freq: BPF cutoff frequencies
% nSamp_predict: number of samples to be predicted (further after the segment lost by filter delay)
% 
% Written by Seonghun Park (s.park7532@gmail.com) according to the article below:
% Mansouri, F., Dunlop, K., Giacobbe, P., Downar, J., & Zariffa, J. (2017). A fast EEG forecasting algorithm for phase-locked transcranial electrical stimulation of the human brain. Frontiers in neuroscience, 11, 401.
%
% Every parameter is set up as in the reference


data_len = length(data)/fs;

ord = 10;  
passband_ripple = 0.5;  
stopband_ripple = 40;   
[A, B, C, D] = ellip(ord, passband_ripple, stopband_ripple, cutoff_freq/(fs/2));
sos = ss2sos(A, B, C, D);
data_filtered = sosfilt(sos, data);

NFFT_padded = 20*fs;

FFT = fft(data_filtered, NFFT_padded)/length(data_filtered);
FFT_amp = 2*abs(FFT);
f = linspace(0, fs/2, NFFT_padded/2+1);


f_resol = fs/NFFT_padded;
n_step = 1/f_resol;

[FFT_amp_target_f, max_idx] = max(FFT_amp(cutoff_freq(1)*n_step+1:cutoff_freq(2)*n_step));
max_freq_idx = max_idx + cutoff_freq(1)*n_step;  % index from full frequency range


[B, A] = ellip(ord, passband_ripple, stopband_ripple, cutoff_freq/(fs/2));
[~, phase_resp] = freqz(B, A, 10*fs, fs);


FFT_phase = angle(FFT);
FFT_phase_target_f = FFT_phase(max_freq_idx);
FFT_phase_target_f_corrected = FFT_phase_target_f + phase_resp(max_freq_idx);


t_prediction = 1/fs:1/fs:data_len+nSamp_predict/fs;
target_f = f(max_freq_idx);
output = FFT_amp_target_f*cos(2*pi*target_f*t_prediction + FFT_phase_target_f_corrected - phase_resp(max_freq_idx));


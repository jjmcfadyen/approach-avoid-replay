function [freq,amp,phase] = getPhase(t,s)
% t = time vector
% s = signal vector

Ts = mean(diff(t));                                         % Sampling Time
Fs = 1/Ts;                                                  % Sampling Frequency
Fn = Fs/2;                                                  % Nyquist Frequency

L  = length(s);
fts = fft(s)/L;                                             % Normalised Fourier Transform
Fv = linspace(0, 1, fix(L/2)+1)*Fn;                         % Frequency Vector
Iv = 1:length(Fv);                                          % Index Vector

amp_fts = abs(fts(Iv))*2;                                   % Spectrum Amplitude
phs_fts = angle(fts(Iv));                                   % Spectrum Phase

% ---
freq = Fv(Fv>0);
amp = amp_fts(Fv>0);
phase = phs_fts(Fv>0);

% %
% 
% subplot(1,2,1)
% plot(freq,amp);
% title(['Max freq = ' num2str(round(freq(amp==max(amp)))) ' Hz'])
% 
% subplot(1,2,2)
% plot(freq,phase);

end
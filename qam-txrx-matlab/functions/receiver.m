function data_bits=receiver(signal_rx,Frame_length,s)

% Variable initialization
Fs_rx=s.Fs_rx; %Sampling rate (Msamples/s)
Rs=s.Rs; %Symbol rate (Msymbols/s)

Delay=s.Delay;
Rolloff=s.Rolloff;

SampSymb_rx=Fs_rx/Rs;
Num_pilot_symbols=s.Num_data_pilots;
Num_data_symbols=s.Num_data_symbols;
pilot_symbols_tx=s.pilot_symbols_tx;

M=s.M;

norm_factor_qam=[1 2 NaN 10 NaN 42];
tx_user=1;

Observables=[];

%% BEGIN: Key symbols
key_bits=user_key(tx_user);
key_symbols=2*key_bits-1; %Modulate BPSK

% Key symbols pulse shaping
Fd = 1;
coef = rcosine(Fd,SampSymb_rx,'sqrt',Rolloff,Delay);
key = rcosflt(key_symbols,Fd,SampSymb_rx,'filter',coef).';
key = key(1:end-SampSymb_rx*Delay+1);
%%END: Key symbols

%% BEGIN: Signal preprocessing

% Suppress DC offset
Long_sec_adc=size(signal_rx,2);
x=signal_rx-kron(mean(signal_rx,2),ones(1,Long_sec_adc));

% Coarse frequency offset estimation
Nfft=2.^(nextpow2(Long_sec_adc) + 4);
% Nfft=2^20;
X = fftshift(fft(x.^2,Nfft));
freqs=((2*pi/Nfft) * [0:(Nfft-1)])*Fs_rx/(2*pi);
mid=ceil(Nfft/2)+1;
freqs(mid:Nfft)=freqs(mid:Nfft)-Fs_rx;		%<-- move [pi,2pi] to [-pi,0]
freqs=fftshift(freqs);
[value,index] = max(abs(X).*(abs(freqs)<=(2*(Fs_rx/SampSymb_rx)*(1 + Rolloff))));
ft = freqs(index)/2;

offset_signal = exp(-1i*2*pi*ft*(0:(Long_sec_adc-1))/Fs_rx);

thr = 0.2e6; %If offset is larger than this do not perform frequency offset correction

% Shift the signal towards f=0
if abs(ft)<thr
    x_freq_corr=x.*offset_signal;
else
    warning off backtrace;
    warning('No frequency correction due to excessive offset');
    warning on backtrace;
    x_freq_corr=x;
end

% Normalization
mx = max(max([abs(real(x_freq_corr)), abs(imag(x_freq_corr))]));
x_freq_corr = x_freq_corr/mx;


% figure;
% plot(freqs,10*log10(abs(X)))
%%END: Signal preprocessing

%% BEGIN: Baseband processing
%%BEGIN: Correlation with key

Rxy = xcorr(x_freq_corr,key);

r=Rxy(Long_sec_adc:end);
r=r(1:end-Frame_length+1);

[mx,n]=max(abs(r));

% figure
% plot(1:length(r),1e-1*abs(r),'r')
% hold on;
% plot(1:length(x_freq_corr),real(x_freq_corr),'b')
% plot(1:length(x_freq_corr),imag(x_freq_corr),'m')
%%END: Correlation with key

%% Extract frame

frame_rx=signal_rx(n:n+Frame_length-1);
frame_rx=frame_rx-mean(frame_rx);

figure(101)
subplot(2,1,1)
plot(real(frame_rx),'b')
xlabel('Time samples')
ylabel('Received signal - real part')
subplot(2,1,2)
plot(imag(frame_rx),'b')
xlabel('Time samples')
ylabel('Received signal - imaginary part')
set(gcf,'WindowStyle','docked')

%% BEGIN: Demodulation

coef=rcosine(Fd,SampSymb_rx,'sqrt',Rolloff,Delay);

% Filter received signal
[frame_rx_filt,aux]=rcosflt(frame_rx,Fd,SampSymb_rx,'Fs/filter',coef);
% Column to vector
frame_rx_filt = frame_rx_filt.';
% Remove filter transients
frame_rx_filt = frame_rx_filt(1+2*Delay*SampSymb_rx:end-2*Delay*SampSymb_rx);

% Get the observations
obs = frame_rx_filt(1:SampSymb_rx:end);
% Suppress the key
obs = obs(1+length(key_symbols):end);

% obs_rot=obs;
[obs_rot,Hg_rot] = estima_frec_MIMO(obs,pilot_symbols_tx,20); %ToDo, rewrite this function

pilot_symbols_rx=obs_rot(1:Num_pilot_symbols);
data_symbols_rx=obs_rot(1+Num_pilot_symbols:end);

Hg = sum(pilot_symbols_rx./pilot_symbols_tx)/Num_pilot_symbols;

% Channel equalization: SISO
pilot_symbols_eq = pilot_symbols_rx/Hg;
data_symbols_eq = data_symbols_rx/Hg;

figure(1001)
plot(pilot_symbols_eq.','r.');hold on
plot(data_symbols_eq.','b.')
%     plot(pilot_symbols_rx.','kx');hold on
%     plot(data_symbols_rx.','mx');
maximo = 1.2*max(max(abs(pilot_symbols_eq)));
axis('square');grid
axis([-maximo maximo -maximo maximo])
grid on;
set(gcf,'WindowStyle','docked')

%Plot eye-diagrams
Num_key_symbols=length(key_symbols);
frame_rx_filt_rot=frame_rx_filt/Hg;
% eyediagram(real(frame_rx_filt_rot(1:SampSymb_rx*Num_key_symbols)),SampSymb_rx,SampSymb_rx,0,'g-')
% eyediagram(real(frame_rx_filt_rot(SampSymb_rx*Num_key_symbols+1:SampSymb_rx*(Num_key_symbols+Num_pilot_symbols))),SampSymb_rx,SampSymb_rx,0,'r-')
% eyediagram(frame_rx_filt_rot(SampSymb_rx*(Num_key_symbols+Num_pilot_symbols)+1:end),SampSymb_rx,SampSymb_rx,0,'b-')

% Decision
QAMdemodulator=modem.qamdemod('M',M,'SymbolOrder','Gray','OutputType','Bit','DecisionType','hard decision');
data_bits=demodulate(QAMdemodulator,sqrt(norm_factor_qam(log2(M)))*data_symbols_eq);
% return

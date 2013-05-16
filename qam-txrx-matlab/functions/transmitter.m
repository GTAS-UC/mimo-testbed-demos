function [frame pilot_symbols data_symbols]=transmitter(data_bits,s)

% Variable initialization
Fs_tx=s.Fs_tx; %Sampling rate (Msamples/s)
Rs=s.Rs; %Symbol rate (Msymbols/s)

Delay=s.Delay;
Rolloff=s.Rolloff;

SampSymb_tx=Fs_tx/Rs;
Num_pilot_symbols=s.Num_data_pilots;
Num_data_symbols=s.Num_data_symbols;

M=s.M;

norm_factor_qam=[1 2 NaN 10 NaN 42];
tx_user=1;

%% BEGIN: Key symbols

key_bits=user_key(tx_user);
key_symbols=2*key_bits-1; %Modulate BPSK

%%END: Key symbols

%% BEGIN: Pilot symbols

Num_pilot_bits=Num_pilot_symbols; %BPSK pilot symbols

pilot_symbols = sign(randn(1,Num_pilot_bits)); %BPSK random pilot bits

%%END: Pilot symbols

%% BEGIN: Data Symbols

QAMmodulator=modem.qammod('M',M,'SymbolOrder','Gray','InputType','Bit');
data_symbols=1/sqrt(norm_factor_qam(log2(M)))*modulate(QAMmodulator,data_bits);

%% Perform pulse shaping
Fd=1;

coef = rcosine(Fd,SampSymb_tx,'sqrt',Rolloff,Delay);

tmp = rcosflt([key_symbols pilot_symbols data_symbols],Fd,SampSymb_tx,'filter',coef).';
frame= tmp/max(abs([real(tmp) imag(tmp)]));

% plot(real(frame),'r')
% hold on;
% plot(imag(frame),'b')

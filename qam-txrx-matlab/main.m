clc; close all; clear;

addpath wsLyrtechGTAS;
addpath functions;
%% Create objects from web service
tx = wsLyrtechGTAS('http://www.gtas.dicom.unican.es/testbed/node1/wsLyrtech.asmx');
rx = wsLyrtechGTAS('http://www.gtas.dicom.unican.es/testbed/node2/wsLyrtech.asmx');

%% Common configuration parameters
f_RF = 5600;

%% Auxiliary functions
ant2chan=@(x)reshape([2*x-1; 2*x],[],1)';
chan2ant=@(x)unique(fix((x-1)/2)+1);

%% TX Configuration parameters

dac_v2 = ~eval(dac_IsVirtex4(tx));

if dac_v2
    bit_file = 'DAC_VirtexII_8_channels_RFFE.bit';interpolation_tx = 1;
else
    bit_file = 'DAC_Virtex4_8_channels_RFFE.bit';interpolation_tx = 1;
%     bit_file = 'DAC_Virtex4_8_channels_RFFE_interpolate_2.bit';interpolation = 2;
    %bit_file = 'DAC_Virtex4_8_channels_RFFE_interpolate_4.bit';interpolation = 4;
end
programmed_dac = dac_FPGAConf(tx,bit_file);

programmed_adc = adc_FPGAConf(tx,'ADC_Virtex4_8_channels_RFFE.bit'); % The ADC-FPGA is programmed since the RFFE is controlled through it

tx_antennas = [2];
tx_channels = ant2chan(tx_antennas);

sampling_freq_TX = 52;
sampling_freq_TX_DAC = sampling_freq_TX*interpolation_tx;
DAC_gain = 255*dac_v2+15*(1-dac_v2);   % 0 to 255 for VirtexII and 0 to 15 for Virtex4
channel_gains_set = dac_ChannelGains(tx,DAC_gain*ones(1,8));

freq_set = dac_SamplingFreq(tx,sampling_freq_TX_DAC);

fc = 24; %Corner freq in MHz
rffe_initialized = rffe_Init(tx,fc,f_RF);

%Switch on needed transceivers
for tr = tx_antennas
    pa_gain_set = rffe_TxGain(tx,20,tr);
end

%Switch off not needed transceivers
for tr = setxor(1:4,tx_antennas)
    pa_gain_unset = rffe_TxGain(tx,-1,tr);
end

%% RX Configuration parameters
bit_file = 'ADC_Virtex4_8_channels_RFFE.bit';decimate = 1;
% bit_file = 'ADC_Virtex4_8_channels_RFFE_decimate_2.bit';decimate = 2;
%bit_file = 'ADC_Virtex4_8_channels_RFFE_decimate_4.bit';decimate = 4;

% FPGA programming
programmed_adc = adc_FPGAConf(rx,bit_file);

rx_antennas = [2];
rx_channels = ant2chan(rx_antennas);

sampling_freq_RX = 52;
sampling_freq_RX_DAC = sampling_freq_RX/decimate;

freq_set = adc_SamplingFreq(rx,sampling_freq_RX_DAC);

fc = 20; %Corner freq in MHz
rffe_initialized = rffe_Init(rx,fc,f_RF);

%Switch on needed transceivers
for tr = rx_antennas
    lna_vga_gain_set = rffe_RxGain(rx,3,0,tr);
end
%Switch off not needed transceivers
for tr = setxor(1:4,rx_antennas)
    lna_vga_gain_gain_unset=rffe_RxGain(rx,-1,-1,tr);
end

%ToDo add external synchronization (move to a function and fix rffe_LOCKED)
adc_WriteFPGAReg(rx,hex2dec('123C'),hex2dec('2'));
%% Actual transmission and reception
%%Parameters

s.Num_data_pilots=256;
s.Num_data_symbols=512;
s.M=4; % M-QAM modulation

s.Fs_tx=52e6;
s.Fs_rx=52e6; %Sampling rate (Msamples/s)
s.Rs=2e6; %Symbol rate (Msymbols/s)

s.Delay=10;
s.Rolloff=0.4;

data_bits=randi([0 1],log2(s.M),s.Num_data_symbols);

%% Transmit
[frame_tx pilot_symbols_tx data_symbols_tx]=transmitter(data_bits,s);

% Save new parameters
s.pilot_symbols_tx=pilot_symbols_tx;
s.data_symbols_tx=data_symbols_tx;

v(1,:) = real(frame_tx);
v(2,:) = imag(frame_tx);

nSamples_tx=length(frame_tx);

if dac_v2
    v_DAC = v*(2^(11)-1);
else
    v_DAC = -v*(2^(11)-1);
end
num_samples = dac_SDRAMWrite(tx,v_DAC); %v_DAC must have one channel per row
channels_configured = dac_ChannelConf(tx,sum(2.^(tx_channels-1)),length(v_DAC));
playback_started = dac_Play(tx);

%% Receive
close all
num_samples = 100000;

channels_configured = adc_ChannelConf(rx,sum(2.^(rx_channels-1)),num_samples);
data = adc_Acquire_SDRAMRead(rx,num_samples);

% Virtual transmission
% data=repmat(v,1,ceil(num_samples/nSamples_tx));
% data=data(:,1:num_samples);

signal_rx=data(1,:)+1j*data(2,:);
data_bits_rx=receiver(signal_rx,nSamples_tx,s);
errors=abs(data_bits_rx-data_bits);
BER=sum(errors(:))/numel(errors(:))

return;
%%
playback_stopped = dac_ChannelDeconf(tx);

for tr = 1:4
    pa_gain_unset = rffe_TxGain(tx,-1,tr);
end
node_Unlock(tx);
node_Unlock(rx);
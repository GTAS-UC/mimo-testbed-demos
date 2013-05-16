clc; close all; clear;

%% Create object from web service
node = 'node1';
tx = wsLyrtechGTAS(['http://www.gtas.dicom.unican.es/testbed/' node '/wsLyrtech.asmx']);

dac_v2 = ~eval(dac_IsVirtex4(tx));
%% Configuration parameters
if dac_v2
    bit_file = 'DAC_VirtexII_8_channels_RFFE.bit';interpolation = 1;
else
    bit_file = 'DAC_Virtex4_8_channels_RFFE.bit';interpolation = 1;    
    %bit_file = 'DAC_Virtex4_8_channels_RFFE_interpolate_2.bit';interpolation = 2;
    %bit_file = 'DAC_Virtex4_8_channels_RFFE_interpolate_4.bit';interpolation = 4;
end

bit_file = 'ADC_Virtex4_8_channels_RFFE.bit';decimate = 1;
%bit_file = 'ADC_Virtex4_8_channels_RFFE_decimate_2.bit';decimate = 2;
%bit_file = 'ADC_Virtex4_8_channels_RFFE_decimate_4.bit';decimate = 4;

f_RF = 5600;

channels = [1 2];
n_channels = length(channels);
transceivers = unique(fix((channels-1)/2)+1);

sampling_freq_TX = 52;
sampling_freq_TX_DAC = sampling_freq_TX*interpolation;
DAC_gain = 255*dac_v2+15*(1-dac_v2);   % 0 to 255 for VirtexII and 0 to 15 for Virtex4

%% Transmitted Signals

w_0 = f/sampling_freq_TX;

[N,D] = rat(1/w_0);

num_samples = N*D*1000;

t = 2*pi*f/sampling_freq_TX*(0:num_samples-1);
t2 = 2*t; 
v1 = sin(t);
v2 = cos(t);
v3 = square(t);
v4 = sawtooth(t);
v5 = diric(t,5);
v6 = (-1)*sawtooth(t);
v7 = sin(t2);
v8 = cos(t2);
v = [v1 ; v2 ; v3 ; v4 ; v5 ; v6 ; v7 ; v8];
v = v(channels,:);
if dac_v2
    v_DAC = v*(2^(11)-1);
else
    v_DAC = -v*(2^(11)-1);
end

% FPGA programming 
programmed_dac = dac_FPGAConf(tx,bit_file)
programmed_adc = adc_FPGAConf(tx,'ADC_Virtex4_8_channels_RFFE.bit') % The ADC-FPGA is programmed since the RFFE is controlled through it

fc = 24; %Corner freq in MHz
rffe_initialized = rffe_Init(tx,fc,f_RF)

%Switch on needed transceivers
for tr = transceivers
    disp(sprintf('Transceiver %d',tr))
    pa_gain_set = rffe_TxGain(tx,30,tr)
end

%Switch off not needed transceivers
for tr = setxor(1:4,transceivers)
    disp(sprintf('Transceiver %d',tr))
    pa_gain_unset = rffe_TxGain(tx,-1,tr)
end


freq_set = dac_SamplingFreq(tx,sampling_freq_TX_DAC)
num_samples = dac_SDRAMWrite(tx,v_DAC) %v_DAC must have one channel per row
channels_configured = dac_ChannelConf(tx,sum(2.^(channels-1)),length(v_DAC))
channel_gains_set = dac_ChannelGains(tx,DAC_gain*ones(1,8))

playback_started = dac_Play(tx)

return;

pause(10)
playback_stopped = dac_ChannelDeconf(tx)

for tr = 1:4
    disp(sprintf('Transceiver %d',tr))
    pa_gain_unset = rffe_TxGain(tx,-1,tr)
end

node_Unlock(tx)
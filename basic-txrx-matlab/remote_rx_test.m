clc; close all; clear;

%% Create object from web service
node = 'node2'
rx = wsLyrtechGTAS(['http://www.gtas.dicom.unican.es/testbed/' node '/wsLyrtech.asmx']);

%% Configuration parameters

bit_file = 'ADC_Virtex4_8_channels_RFFE.bit';decimate = 1;
%bit_file = 'ADC_Virtex4_8_channels_RFFE_decimate_2.bit';decimate = 2;
%bit_file = 'ADC_Virtex4_8_channels_RFFE_decimate_4.bit';decimate = 4;

f_RF = 5600;

channels = [1 2];
n_channels = length(channels);
transceivers = unique(fix((channels-1)/2)+1);

num_samples = 1000;

sampling_freq_RX = 52;
sampling_freq_RX_DAC = sampling_freq_RX/decimate;

% FPGA programming
dac_v2 = ~eval(dac_IsVirtex4(rx));
if dac_v2
    programmed_dac = dac_FPGAConf(rx,'DAC_VirtexII_8_channels_RFFE.bit')
else
    programmed_dac = dac_FPGAConf(rx,'DAC_Virtex4_8_channels_RFFE.bit')
end
programmed_adc = adc_FPGAConf(rx,bit_file)

fc = 20 %Corner freq in MHz
rffe_initialized = rffe_Init(rx,fc,f_RF)

%Switch on needed transceivers
for tr = transceivers
    disp(sprintf('Transceiver %d',tr))
    lna_vga_gain_set = rffe_RxGain(rx,3,0,tr)
end
%Switch off not needed transceivers
for tr = setxor(1:4,transceivers)
    disp(sprintf('Transceiver %d',tr))
    lna_vga_gain_gain_unset=rffe_RxGain(rx,-1,-1,tr)
end

freq_set = adc_SamplingFreq(rx,sampling_freq_RX_DAC)
channels_configured = adc_ChannelConf(rx,sum(2.^(channels-1)),num_samples)
data = adc_Acquire_SDRAMRead(rx,num_samples);

plot(data.')

node_Unlock(rx)
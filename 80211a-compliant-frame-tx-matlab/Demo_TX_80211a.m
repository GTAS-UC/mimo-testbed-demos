clc; close all; clear;

%% Configure transmitter

%Parameters
node = 'node1';
f_RF = 5600; % 5.6 GHz = Channel 120
channels = [1 2]; %DAC channels

tx = wsLyrtechGTAS(['http://gtas.unican.es/testbed/' node '/wsLyrtech.asmx']);

switch lower(node)
    case 'node1'
%         bit_file = 'DAC_VirtexII_8_channels_RFFE.bit';interpolation = 1;
        bit_file = 'DAC_VirtexII_8_channels_RFFE_interpolate_2_fir.bit';interpolation = 2;
        DAC_gain=255;
        sampling_freq_TX_DAC=40;
    case {'node2','node3'}
        bit_file = 'DAC_Virtex4_8_channels_RFFE.bit'; interpolation = 1;
%         bit_file = 'DAC_Virtex4_8_channels_RFFE_interpolate_2.bit'; interpolation = 2;
%         bit_file = 'DAC_Virtex4_8_channels_RFFE_interpolate_4.bit'; interpolation = 4;
        DAC_gain=15;
        sampling_freq_TX_DAC=52;
    otherwise
        error('Wrong node name')
end

n_channels = length(channels);
transceivers = unique(fix((channels-1)/2)+1);

programmed_dac = dac_FPGAConf(tx,bit_file)
programmed_adc = adc_FPGAConf(tx,'ADC_Virtex4_8_channels_RFFE.bit') % The ADC-FPGA is programmed since the RFFE is controlled through it

fc = 18; %Corner freq in MHz
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

freq_set = eval(dac_SamplingFreq(tx,sampling_freq_TX_DAC))

%% Generate frame
% Frame parameters
RATE=6;

mac_framectrl = '0802';          % this is a data frame
mac_duration  = '002e';          % is calculated further down
mac_address1  = 'ffffffffffff';  % destination mac-address
mac_address2  = '001122334455';  % BSSID mac-address
mac_address3  = '001122334455';  % source mac-address
mac_seqctrl   = '0000';          % 0 sequence and 0 fragment number
mac_crc       = '00000000';      % is calculated further down
mac_header = strcat(mac_framectrl,mac_duration,mac_address1,mac_address2,mac_address3,mac_seqctrl);

% LONG=50;
% msg=randint(1,LONG,[0 15]); %random message
msg='Esto es una trama 802.11a'; %user-defined message
msg_hex = reshape(dec2hex(msg).',1,[]); %ASCII message to HEX message
crc = crc32([mac_header msg_hex]); %Compute valid CRCs
% crc='00000000'; %Force a wrong CRC

% Build data frame
data_frame_pre_interp = data_frame_80211a(RATE,[mac_header msg_hex crc]);

% Interpolate, if needed
[p q]=rat(freq_set/20/interpolation);
data_frame_interp=resample(data_frame_pre_interp,p,q);
% data_frame_interp=data_frame_pre_interp;

v=[real(data_frame_interp); imag(data_frame_interp)];

v_DAC = v*(2^(11)-1);

%%% Upload signal
num_samples = dac_SDRAMWrite(tx,v_DAC) %v_DAC must have one channel per row
channels_configured = dac_ChannelConf(tx,sum(2.^(channels-1)),length(v_DAC))
channel_gains_set = dac_ChannelGains(tx,DAC_gain*ones(1,8))

%% Enable single-shot playback
singleshot = dac_PlaybackMode(tx,'single-shot')

return
%% Enable continuous playback
continuous = dac_PlaybackMode(tx,'continuous')

%% Transmit now!
playback_started = dac_Play(tx)

%%
channels_deconfigured = dac_ChannelDeconf(tx)

%% Configure laptop to receive frames
%{
sudo service network-manager stop
sudo airmon-ng start wlan0 120
sudo wireshark
Analyze/Enabled Protocols.../802.11 Radiotap & IEEE 802.11
%}
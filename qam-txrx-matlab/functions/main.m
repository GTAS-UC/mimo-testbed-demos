clc;
clear;
close all;

[frame_tx pilot_symbols_tx data_symbols_tx]=transmitter();


signal_rx=circshift(repmat(frame_tx,1,10),[0 randi([1 length(frame_tx)])]);
noise=randn(size(signal_rx))+1j*randn(size(signal_rx));
signal_rx=signal_rx+0.05*noise;

h=randn(1)+1j*randn(1);
signal_rx=signal_rx*h;
% plot(real(signal_rx),'r')
% hold on;
% plot(imag(signal_rx),'b')


figure(101)
plot(real(frame_tx))

s.pilot_symbols_tx=pilot_symbols_tx;
s.data_symbols_tx=data_symbols_tx;

Frame_length=length(frame_tx);
receiver(signal_rx,Frame_length,s);


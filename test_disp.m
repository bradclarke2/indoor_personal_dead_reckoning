clear all
clc

dt = 0.01; % seconds per sample 
N = 512; % number of samples 
t = 0 : dt : (N-1)*dt; % in seconds 
wave_freq = 1; % in Hertz 
acc_time_data = sin(2*pi*wave_freq*t);

data = acc2disp(acc_time_data, 0.01);
function[Filt] = butterworthLowV1(datatext, nocol, Fs, Fc)
%%%%%%%%%%%%%%%%%%%%%STEP TWO%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DG: Changed to take into account normalised cutoff freq (2*Fc/Fs)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        1D- Fast Fourier Transform (FFT)
% Inputs: TempData.Raw
% Outputs: Data.LF, Data.HF
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% -----------------------------------------------------------------------
%                 DEFINE DATA FILE LOCATION AND FILE NAME
% -----------------------------------------------------------------------
%tic

% Variables.FileName = 'datatext.txt';
% Variables.Path = 'C:\Documents and Settings\mmses2\My Documents\Back up Sept 09\PhD Back Up\14th july testing\14th July testing_Node data\FFT_ButterFilter';
% cd(Variables.Path);
% 

%datatext=converted data values, i.e. output from conversionfouraxes m file

Data2=datatext;

%Prompt input for filter order & cut-off frequency
% evalResponse1 = input('please input the required filter order :');
% evalResponse2 = input('please input a filter frequency :');

%sampling rate=0.04s
timet(1)=0.02;
for aa=2:length(Data2)
    timet(aa)=timet(aa-1)+0.02;
end

%u=number of columns that you are filtering
for u=1:nocol;
Data = datatext(:,u);
%Data = jamesDone;
Datatwo=datatext;

% -----------------------------------------------------------------------
%                     DEFINE PROGRAM CONSTANTS
% -----------------------------------------------------------------------
Variables.ZeroPadSize = 180000;          % Zero pad data to defined size
Variables.SampleFreq = Fs;              % Sample Frequency


% -----------------------------------------------------------------------
%                     DEFINE PROGRAM VARIABLES
% -----------------------------------------------------------------------

Variables.n = 5;                        % Filter order
Variables.FilterFreq = Fc;               %Filter frequency for dives
%Variables.FilterFreq = 6;              % Filter Frequency for breaststroke
%Variables.FilterFreq = 2;               % Filter Frequency for front crawl, backstroke
%Variables.FilterFreq = 8;              % Filter frequency for butterfly
%Variables.FilterFreq = 1;              % Lap count for butterfly, backstroke
%Variables.FilterFreq = 0.6;            % Lap count 


% -----------------------------------------------------------------------
%                 GENERATE 1D LOW PASS BUTTERWORTH FILTER
% -----------------------------------------------------------------------

% Generate Frequency Axis to plot results against
if mod(Variables.ZeroPadSize,2)==0
    k=-Variables.ZeroPadSize/2:Variables.ZeroPadSize/2-1; % N even
else
    k=-(Variables.ZeroPadSize-1)/2:(Variables.ZeroPadSize-1)/2; % N odd
end
T=Variables.ZeroPadSize/Variables.SampleFreq;
freq=(2*k)/T;


% figure(1)
% 
% %Figure showing effect of changing filter order on filter response
% subplot(3,2,1)
% for n=1:1:10;
% [B,A] = butter(n,(Variables.FilterFreq/(Variables.SampleFreq)),'Low');
% Butterworth=abs(freqz(B,A,Variables.ZeroPadSize/2));
% Mask.Butterworth=[(flipdim(Butterworth,1)); Butterworth];
% plot(freq,Mask.Butterworth);
% title('Butterworth Filter for Range of Filter Orders');
% xlabel('Frequency (Hz)');
% ylabel('Amplitude');
% hold on
% end
% % hold off



% Determine Butterworth filter coefficients using program variables
[B,A] = butter(Variables.n,(2*Variables.FilterFreq/(Variables.SampleFreq)),'Low');
b=B;
a=A;
Butterworth=abs(freqz(B,A,Variables.ZeroPadSize/2));
Mask.Butterworth=[(flipdim(Butterworth,1)); Butterworth];

% % Figure showing filter used in analysis
% subplot(3,2,2)
% plot(freq,Mask.Butterworth);
% title('Filter Applied to Data');
% xlabel('Frequency (Hz)');
% ylabel('Amplitude');
% hold on




% -----------------------------------------------------------------------
%                  COMPLETE FAST FOURIER TRANSFORM ON INPUT DATA
% -----------------------------------------------------------------------

% Check for NaNs in raw data and set to zero
TempData.Raw=Data;
TempData.Raw(find(isnan(TempData.Raw)==1))=0;

% Complete the FFT on input data, zero-padding to nearest base2 number
TempData.FFT = fftshift(fft(TempData.Raw,Variables.ZeroPadSize));


% Apply Mask to obtain low frequency component
TempData.LF = Mask.Butterworth.*TempData.FFT;

% Transfer from Fourier Domain back into Time Domian
% Take only real component of data following return to time domain
TempData.LF = real(ifft(ifftshift(TempData.LF)));


% Remove ZeroPadded Data and save results
OutputData.Raw = TempData.Raw(1:length(TempData.Raw));
OutputData.LF = TempData.LF(1:length(TempData.Raw));
OutputData.HF = OutputData.Raw-OutputData.LF;

% Store complete zero-padded FFT data for further analysis
OutputData.FFT = TempData.FFT(1:length(TempData.FFT));

Filt(:,u) = OutputData.LF;
% if u==1
%     Filt1 = OutputData.LF;
% end
% if u==2
%     Filt2 = OutputData.LF;
% end
% if u==3
%     Filt3 = OutputData.LF;
% end
% if u==4
%     Filt4 = OutputData.LF;
% end
% if u==5
%     Filt5 = OutputData.LF;
% end
% if u==6
%     Filt6 = OutputData.LF;
% end
u=u+1;
end
%toc

% Clear used variables
clear TempData Butterworth
clear A B T  k n


% % Plot Frequency Spectrum
% % subplot(3,2,3:4)
% % stem(freq,abs(OutputData.FFT), '.');
% % title('Processed Data Frequency Spectrum');
% % xlabel('Frequency (Hz)');
% % ylabel('Amplitude');


% % Plot Filtered Results
% subplot
% % plot(OutputData.LF,'-r')
% % hold on
% plot(OutputData.HF,'-b')
% hold on
% plot(OutputData.Raw,'-k')
% title('Filter Applied to Data Showing Raw, High and Low Frequencies')
% ylabel('Amplitude');
% % hold on
% % plot(real(ifft(lowPass)),'m');
% hold off
% 
% figure(10)
% plot(OutputData.LF,'-r')
% hold on
% plot(OutputData.Raw,'-k')
% hold on
% plot(real(ifft(lowPass)),'m');


% figure(11)
% plot(OutputData.LF,'-r')
% hold on
% plot(OutputData.HF,'-b')
% hold on
% plot(OutputData.Raw,'-k')
% title('Filter Applied to Data Showing Raw, High and Low Frequencies')
% ylabel('Amplitude');
% 
% figure(12)
% plot(timet, FilteredData)
% %title('Filtered applied to data, showing filtered values', 'fontsize', 16)
% ylabel('acceleration (g)', 'fontsize', 16)
% xlabel('time (s)','fontsize', 16)
% 
% figure(13)
% plot(timet,Datatwo)
% title('Raw data','fontsize', 16)
% ylabel('acceleration (g)','fontsize', 16)
% xlabel('time (s)','fontsize', 16)


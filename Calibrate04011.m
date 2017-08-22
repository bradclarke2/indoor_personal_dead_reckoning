function [calD] = Calibrate04011(datainput);

%convert raw accelerometer values to m/s/s
aax = datainput(:,1);
aay = datainput(:,2);
aaz = datainput(:,3);
g1 = datainput(:,4);
g2 = datainput(:,5);
g3 = datainput(:,6);
% Node used : 04011
%Calibration data on Excel sheet 05/03/2013
%Accelerometers
%Column1
calD(:,1)=(aax-521.7832646667)/100.5000450000;
%Column2
calD(:,2)=(aay-508.1082205000)/101.6309750000;
%Column3
calD(:,3)=(aaz-510.9889905000)/103.4187365000;
%column4
gx1=((g1-492.1892776715)/15.97735815);
gxx = mode(gx1(1:10));%ensure gyro is centred around zero by subtracting the mode
calD(:,4)=gx1-gxx;
%column5
gy1=((g2-449.8833732043)/14.9881465);
gyy = mode(gy1(1:10));
calD(:,5)=gy1-gyy;
%column6
gz1 = ((g3-504.3627774)/20.21534712);%adjusted again as freq was 25hz for angle intgrateion 38.862);%re-adjusted as per clean 1 Dgordon 11.04.2011 (was 15.212);
gzz = mode(gz1(1:10));
calD(:,6) = gz1-gzz;
end

function [calD] = Calibrate04008(datainput);

%convert raw accelerometer values to m/s/s
aax = datainput(:,1);
aay = datainput(:,2);
aaz = datainput(:,3);
g1 = datainput(:,4);
g2 = datainput(:,5);
g3 = datainput(:,6);
% Node used : 04008
%Calibration data on Excel sheet 05/03/2013
%Accelerometers
%Column1
calD(:,1)=(aax-517.1956743333)/101.9314276667;
%Column2
calD(:,2)=(aay-495.7909146667)/103.4164306667;
%Column3
calD(:,3)=(aaz-502.9201663333)/102.9039090000;
%column4
gx1=((g1-491.31793617)/16.35043968);
gxx = mode(gx1(1:10));%ensure gyro is centred around zero by subtracting the mode
calD(:,4)=gx1-gxx;
%column5
gy1=((g2-454.9331892787)/14.9757539);
gyy = mode(gy1(1:10));
calD(:,5)=gy1-gyy;
%column6
gz1 = ((g3-506.0676096739)/20.40949625);%adjusted again as freq was 25hz for angle intgrateion 38.862);%re-adjusted as per clean 1 Dgordon 11.04.2011 (was 15.212);
gzz = mode(gz1(1:10));
calD(:,6) = gz1-gzz;
end
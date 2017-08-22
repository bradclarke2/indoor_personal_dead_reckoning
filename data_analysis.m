%% Clearing of previous data

clear all
clc

%% Selection of .csv file and the useful data
 
%files = dir('*.csv');
%for file = files'
    [FileName,PathName] = uigetfile({'*'},'Load Turn File');
    %get filename and use this to append to graphs and figures
    %FileName = pwd;
    %FileName = strcat(FileName,'\',file.name);
    [pathstr, name, ext] = fileparts(FileName); % inputdlg('What is the name of this test?','Test Name Input');
    delimiter = ',';
    FileName = fullfile(PathName,[name ext]);
    %Format string for each line of text:
    formatSpec = '%f%f%f%f%f%f%f%*s%*s%[^\n\r]';
    % Open the text file.
    fileID = fopen(FileName,'r');
    % Read columns of data according to format string.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    % Close the text file.
    fclose(fileID);
    %remove unnecesary columns from raw CSV file
    raw_data = [dataArray{:, 2} dataArray{:, 3} dataArray{:, 4} dataArray{:, 5} dataArray{:, 6} dataArray{:, 7}];
    %Calibrate data using 04011 node calibration file
    calibrated_data = Calibrate04011(raw_data);
    
%% Calibration of data

    %seperate accel and gyro
    raw_accel = [raw_data(:, 1) raw_data(:, 2) raw_data(:, 3)];
    raw_gyro = [raw_data(:, 4) raw_data(:, 5) raw_data(:, 6)];
    cal_accel = [calibrated_data(:, 1) calibrated_data(:, 2) calibrated_data(:, 3)];
    cal_gyro = [calibrated_data(:, 4) calibrated_data(:, 5) calibrated_data(:, 6)];

%% Sizing of Graphs       
    max_y_axis = 25;
    max_x_axis = size(cal_accel,1);
    max_gyro_axis = 3;
    

%%  Calibrated accel and gyro graphs

    figure('name',strcat(name,' - Calibrated 3 axis acceleration'));
    plot(cal_accel);
    xlabel('Frame (delta t = 0.02s)'); % x-axis label
    ylabel('g'); % y-axis label
    axis([0 ,max_x_axis,-1.5,1.5]);
    title (strcat(name,' - Calibrated 3 axis acceleration'));
    legend('X Axis', 'Y Axis', 'Z Axis');

    figure('name',strcat(name,' - Calibrated 3 axis gyro'));
    plot(cal_gyro);
    xlabel('Frame (delta t = 0.02s)'); % x-axis label
    ylabel('Rad/s'); % y-axis label
    axis([0 ,max_x_axis,-max_gyro_axis,max_gyro_axis]);
    title (strcat(name,' - Calibrated 3 axis gyro'));
    legend('X Axis', 'Y Axis', 'Z Axis');

    %% Calculation producing the normalised accelerations

    normx = (cal_accel(:,1)+1)*9.81;
    normy =cal_accel(:,2)*9.81;
    normz = cal_accel(:,3)*9.81;
    gyrox = cal_gyro(:,1);
    gyroy = cal_gyro(:,2);
    gyroz = cal_gyro(:,3);
    
    %% Frequency and time step calculations

    Fs = 50;
    T = 1/Fs;
    
    %% Fast fourier transform
    L = max_x_axis;
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    Y = fft(normx,NFFT)/L;
    f = Fs/2*linspace(0,1,NFFT/2+1);
    
    

    figure;
    % Plot single-sided amplitude spectrum.
    plot(f,2*abs(Y(1:NFFT/2+1))) 
    title('Single-Sided Amplitude Spectrum of y(t)')
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')
    
    %% Filter
    filterx = butterworthLowV1(normx, 1, 50, 15);
    filtery = butterworthLowV1(normy, 1,50, 15);
    
       
        %%displacement = acc2disp(norm,0.02);
    figure('name',strcat(name,' - Filtered  X axis acceleration'),'units','normalized','outerposition',[0 0 1 1]);
    plot(filterx);
    xlabel('Frame (delta t = 0.02s)'); % x-axis label
    ylabel('m/s2'); % y-axis label
    axis([0 ,max_x_axis,-max_y_axis,max_y_axis]);
    title (strcat(name,' - Filtered, Calibrated & Normalised X axis acceleration'));
    
%     figure('name',strcat(name,' - Calibrated & Normalised  X axis acceleration'),'units','normalized','outerposition',[0 0 1 1]);
%     plot(normx);
%     xlabel('Frame (delta t = 0.02s)'); % x-axis label
%     ylabel('m/s2'); % y-axis label
%     axis([0 ,max_x_axis,-max_y_axis,max_y_axis]);
%     title (strcat(name,' - Calibrated & Normalised X axis acceleration'));
% 
%     norm_accel = [normx normy normz];
% 
%     figure('name',strcat(name,' - Calibrated & Normalised(X) 3 axis acceleration'));
%     plot(norm_accel);
%     xlabel('Frame (delta t = 0.02s)'); % x-axis label
%     ylabel('m/s2'); % y-axis label
%     axis([0 ,max_x_axis,-max_y_axis,max_y_axis]);
%     title (strcat(name,' - Calibrated & Normalised(X) 3 axis acceleration'));
%     legend('X Axis', 'Y Axis', 'Z Axis');
% 
%     
%     %%displacement = acc2disp(norm,0.02);
%     figure('name',strcat(name,' - Calibrated & Normalised  X axis acceleration'),'units','normalized','outerposition',[0 0 1 1]);
%     plot(normx);
%     xlabel('Frame (delta t = 0.02s)'); % x-axis label
%     ylabel('m/s2'); % y-axis label
%     axis([0 ,max_x_axis,-max_y_axis,max_y_axis]);
%     title (strcat(name,' - Calibrated & Normalised X axis acceleration'));
% 
% 
    figure('name',strcat(name,' - Calibrated Y axis acceleration'),'units','normalized','outerposition',[0 0 1 1]);
    plot(normy);
    xlabel('Frame (delta t = 0.02s)'); % x-axis label
    ylabel('m/s2'); % y-axis label
    axis([0 ,max_x_axis,-max_y_axis,max_y_axis]);
    title (strcat(name,' - Calibrated Y axis acceleration'));

    figure('name',strcat(name,' - Calibrate Z axis acceleration'),'units','normalized','outerposition',[0 0 1 1]);
    plot(normz);
    xlabel('Frame (delta t = 0.02s)'); % x-axis label
    ylabel('m/s2'); % y-axis label
    axis([0 ,max_x_axis,-max_y_axis,max_y_axis]);
    title (strcat(name,' - Calibrated Z axis acceleration'));
    
    %% Calculating Jerk
    
    for n = 1 : (length(normx)-1)
       jerkx(n) = (normx(n+1)-normx(n))/T;
       jerky(n) = (normy(n+1)-normy(n))/T;
       jerkz(n) = (normz(n+1)-normz(n))/T;
    end
    
    jerkx = jerkx';
    jerky = jerky';
    jerkz = jerkz';
    
    %% Jerk Plots
    
    figure('name',strcat(name,' - Calibrate X axis Jerk'),'units','normalized','outerposition',[0 0 1 1]);
    plot(jerkx);
    xlabel('Frame (delta t = 0.02s)'); % x-axis label
    ylabel('m/s3'); % y-axis label    
    title (strcat(name,' - Calibrated X axis Jerk')); 
    
      figure('name',strcat(name,' - Calibrate Y axis Jerk'),'units','normalized','outerposition',[0 0 1 1]);
    plot(jerky);
    xlabel('Frame (delta t = 0.02s)'); % x-axis label
    ylabel('m/s3'); % y-axis label    
    title (strcat(name,' - Calibrated Y axis Jerk')); 
    
      figure('name',strcat(name,' - Calibrate Z axis Jerk'),'units','normalized','outerposition',[0 0 1 1]);
    plot(jerkz);
    xlabel('Frame (delta t = 0.02s)'); % x-axis label
    ylabel('m/s3'); % y-axis label    
    title (strcat(name,' - Calibrated Z axis Jerk')); 
    
    
%% Gyroscope Plots

    figure('name',strcat(name,' - Calibrated X axis gyro'),'units','normalized','outerposition',[0 0 1 1]);
    plot(gyrox);
    xlabel('Frame (delta t = 0.02s)'); % x-axis label
    ylabel('Rad/s'); % y-axis label
    axis([0 ,max_x_axis,-max_gyro_axis,max_gyro_axis]);
    title (strcat(name,' - Calibrated X axis gyro'));

    figure('name',strcat(name,' - Calibrated Y axis gyro'),'units','normalized','outerposition',[0 0 1 1]);
    plot(gyroy);
    xlabel('Frame (delta t = 0.02s)'); % x-axis label
    ylabel('Rad/s'); % y-axis label
    axis([0 ,max_x_axis,-max_gyro_axis,max_gyro_axis]);
    title (strcat(name,' - Calibrated Y axis gyro'));

    figure('name',strcat(name,' - Calibrated Z axis gyro'),'units','normalized','outerposition',[0 0 1 1]);
    plot(gyroz);
    xlabel('Frame (delta t = 0.02s)'); % x-axis label
    ylabel('Rad/s'); % y-axis label
    axis([0 ,max_x_axis,-max_gyro_axis,max_gyro_axis]);
    title (strcat(name,' - Calibrated Z axis gyro'));
    
    %% Calculating angle
    anglex = 0;
    angley = 0;
    anglez = 0;
    
    for n = 1: (length(gyrox)-1)
       x_angle =((gyrox(n)+ gyrox(n+1))/2)*T;
       anglex(n+1) = anglex(n) + x_angle; 
       
       y_angle =((gyroy(n)+ gyroy(n+1))/2)*T;
       angley(n+1) = angley(n) + y_angle; 
       
       z_angle =((gyroz(n)+ gyroz(n+1))/2)*T;
       anglez(n+1) = anglez(n) + z_angle; 
    end
 
    
    anglex = (anglex') * 180/pi;
    angley = (angley') * 180/pi;
    anglez = (anglez') * 180/pi;

    figure('name',strcat(name,' - Calibrated X axis angular dsiplacement'),'units','normalized','outerposition',[0 0 1 1]);
    plot(anglex);
    xlabel('Frame (delta t = 0.02s)'); % x-axis label
    ylabel('Degrees'); % y-axis label
    title (strcat(name,' - Calibrated X axis angular displacement'));
    
    figure('name',strcat(name,' - Calibrated Y axis angular displacement'),'units','normalized','outerposition',[0 0 1 1]);
    plot(angley);
    xlabel('Frame (delta t = 0.02s)'); % x-axis label
    ylabel('Degrees'); % y-axis label
    title (strcat(name,' - Calibrated Y axis angular displacement'));
    
    figure('name',strcat(name,' - Calibrated Z axis angular displacement'),'units','normalized','outerposition',[0 0 1 1]);
    plot(anglez);
    xlabel('Frame (delta t = 0.02s)'); % x-axis label
    ylabel('Degrees'); % y-axis label
    title (strcat(name,' - Calibrated Z axis angular displacement'));
    
    %% Calculating velocity profile
%     
%     velocity_x = 0;
%     velocity_y = 0;
%     velocity_z = 0;
%     
%     for n = 1:length(normx)
%         xvelo = T
%         
%     end
% 
%     h = get(0,'children');
%     for i=1:length(h)
%       saveas(h(i), get(h(i),'Name'), 'png');
%     end

%     close all
%end




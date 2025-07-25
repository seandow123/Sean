%--------------------------------------------------------------------------
% Project Propulsion extraction
%   The purpose of the function is to calculate the propulsion data using
%   the legsys sensor
% Author:  Mohsen Zahiri
% Code Description:
%   Version 1.0
%   This code will load in the legsys sensor dat and use the algorithm to
%   detect the propulsion parameters. The algorithm uses the accelometer
%   data to find the propulsion points
%--------------------------------------------------------------------------
%% Initialization
clear all; clc; home;
warning off
option_sensor   = 1;          % 1: Use only Wrist Sensor 2: Use Wrist + Upper Arm Sensor
format short g

%% Directory Location
% Code Directory
CodeDir = fullfile(pwd);
% Raw data directory
% RawDataDir  = 'Z:\Projects BCM\H-38994 BLANKET BCM\Studies\Fscan Frailty\Data\Converted\Gait';

RawDataDir  = 'Z:\Projects BCM\H-38994 BLANKET BCM\Studies\Fscan Frailty\Data\Raw Sensor data\Gait';

%% Load Data
% Load Data File
% Location of the data file
cd(RawDataDir);
directory = [pwd,filesep];
% Opend dialog box to select the folder to analyze.
%   Folder could be analyzied individually or by batch (+1)
listDir = dir(['*']);    % List of Patient Data (starting with PAD*)

% Select Data Folder(s) for Analysis
% s: selection - vector of indices of the selected strings (lenght 1
% in single selection mode. Will be empty([]) when OK is 0.  OK is 1 if you
% push the OK button, and 0 if you push the cancel button in the
% figure
[s,ok] = listdlg('PromptString','Select a folder:', 'ListString',{listDir.name});
clear ok

%% Loop for Analyzing Each Data Folder
%   The code is modified to account for follows up inside each other
%   Optiontime (Baseline - etc.)

PatientIndexData = s(1);

%% define the height of the subject

% Since we have not manage all data, for current run, we have to manually import
%     the height of the subject in the pop up window
%     if the subejcts are from the OHI study, we can find their height from 
%         "Z:\Projects BCM\H-38994 BLANKET BCM\Studies\Fscan Frailty\Data\Raw Sensor data\OHi_Points.xlsx"
%         
%         For the control group we have to look at the check list

% OHI study look at Z:Projects BCM:H-38994 BLANKET BCM:Studies:
% Fscan Frailty:Data:Raw Sensor data:OHi_Points.xlsx, 
% Control group look at the checl list ';
    prompt = 'Enter the Height of the subbject (cm): ';

Height = input(prompt)

%%  load the gait parameters


dir_current = fullfile([directory listDir(PatientIndexData).name],filesep);
% clear PatientIndexData
disp([' '])
disp(['----- Analyzing ', dir_current(end-24:end-2),' ------'])

% Intialize
SensorData = {};
[SensorData, samp_rate,~] = LoadLEGSysRawData5Sensors(dir_current);


%% Raw and Filter Gyroscope Data
%{
        For propulsiond data, two sensors were placed on the left and right
        shin.  Considered Right sensor(1) and Left sensor(2) - however
        would need to recheck data to verify this, but this is to keep the
        same consistent data with previous code
        % patient with negative signal (remove negative sign)
        * Whether the sesnor setup is 5 or 2 from the Legsys will be on S1
        and S3
%}
% Check if there is 2 sensor or 5 sensor set up was used
len = length(SensorData(:,12));
Temp1 = length(find(SensorData(:,12) ==0));
Temp2 = length(find(SensorData(:,13) ==0));
Temp3 = length(find(SensorData(:,14) ==0));

if ( (Temp1==Temp2) & (Temp2==Temp3) & (Temp3==len))
    SensorSetup = 2;
    display(['Sensor setup = ',num2str(SensorSetup)])
else
    SensorSetup = 5;
    display(['Sensor setup = ',num2str(SensorSetup)])
end
clear Temp1 Temp2 Temp3 len

%% filtering and resampling
Gyro_test_RS = filtfilt([1,-1],[1,-0.995],-SensorData(:,4));
Gyro_test_LS = filtfilt([1,-1],[1,-0.995],-SensorData(:,24));

Gyro_test_RS_resample = resample(Gyro_test_RS,200,samp_rate);
Gyro_test_LS_resample = resample(Gyro_test_LS,200,samp_rate);
Fs = 200;
filterorder = 7;
filtercutoff = 15/(Fs/2);
filtertype = 'low';
[b,a] = butter(filterorder,filtercutoff,filtertype);
clear Fs filterorder filtercutoff filtertype
% Right shin
Gyro_St_RS_r = filtfilt(b,a,Gyro_test_RS_resample);
Gyro_St_LS_r = filtfilt(b,a,Gyro_test_LS_resample);
Gyro_St_ss = [Gyro_St_RS_r, Gyro_St_LS_r];
res_St     = GaitAnalyze(Gyro_St_ss, Height);

%% Define the direction of the sensors
if(SensorSetup == 2 || SensorSetup ==5)
    % Sensor 1 or "Right shin"
    Gyro_St_RS(:,1) =  SensorData(:,2);	% x for s1
    Gyro_St_RS(:,2) =  SensorData(:,3);	% y for s1
    Gyro_St_RS(:,3) = -SensorData(:,4);	% z for s1
    Acc_St_RS(:,1) =  SensorData(:,5);	% x for s1
    Acc_St_RS(:,2) =  SensorData(:,6);	% y for s1
    Acc_St_RS(:,3) =  SensorData(:,7);	% z for s1
    
    % Sensor 2 or "Left Shin"
    Gyro_St_LS(:,1) =  SensorData(:,22);	% x for s3
    Gyro_St_LS(:,2) =  SensorData(:,23);	% y for s3
    Gyro_St_LS(:,3) = -SensorData(:,24);	% z for s3
    Acc_St_LS(:,1) =  SensorData(:,25);	% x for s1
    Acc_St_LS(:,2) =  SensorData(:,26);	% y for s1
    Acc_St_LS(:,3) =  SensorData(:,27);	% z for s1
else
    error('Wrong sensor setup')
end

% Filter data
% Right shin
Gyro_St_RS_x =  filtfilt([1,-1],[1,-0.995],-SensorData(:,2));
Gyro_St_RS_y =  filtfilt([1,-1],[1,-0.995],-SensorData(:,3));
Gyro_St_RS   =  filtfilt([1,-1],[1,-0.995],-SensorData(:,4));	% z
Acc_St_RS_x =  filtfilt([1,-1],[1,-0.995],-SensorData(:,5));
Acc_St_RS_y =  filtfilt([1,-1],[1,-0.995],-SensorData(:,6));
Acc_St_RS_z   =  filtfilt([1,-1],[1,-0.995],-SensorData(:,7));	% z
% Left shin
Gyro_St_LS_x =  filtfilt([1,-1],[1,-0.995],-SensorData(:,22));
Gyro_St_LS_y =  filtfilt([1,-1],[1,-0.995],-SensorData(:,23));
Gyro_St_LS   =  filtfilt([1,-1],[1,-0.995],-SensorData(:,24));  % z
Acc_St_LS_x =  filtfilt([1,-1],[1,-0.995],-SensorData(:,25));
Acc_St_LS_y =  filtfilt([1,-1],[1,-0.995],-SensorData(:,26));
Acc_St_LS_z   =  filtfilt([1,-1],[1,-0.995],-SensorData(:,27));	% z
clear StartPt EndPt

% Resampling
Gyro_St_RS_r   = resample(Gyro_St_RS,200,samp_rate);
Gyro_St_RS_x_r = resample(Gyro_St_RS_x,200,samp_rate);
Gyro_St_RS_y_r = resample(Gyro_St_RS_y,200,samp_rate);
Acc_St_RS_z_r   = resample(Acc_St_RS_z,200,samp_rate);
Acc_St_RS_x_r = resample(Acc_St_RS_x,200,samp_rate);
Acc_St_RS_y_r = resample(Acc_St_RS_y,200,samp_rate);

Gyro_St_LS_r   = resample(Gyro_St_LS,200,samp_rate);
Gyro_St_LS_x_r = resample(Gyro_St_LS_x,200,samp_rate);
Gyro_St_LS_y_r = resample(Gyro_St_LS_y,200,samp_rate);
Acc_St_LS_z_r   = resample(Acc_St_LS_z,200,samp_rate);
Acc_St_LS_x_r = resample(Acc_St_LS_x,200,samp_rate);
Acc_St_LS_y_r = resample(Acc_St_LS_y,200,samp_rate);
clear Gyro_St_RS Gyro_St_RS_x Gyro_St_RS_y
clear Gyro_St_LS Gyro_St_LS_x Gyro_St_LS_y

% Low Pass filtering
% Filter parameters
Fs = 200;
filterorder = 7;
filtercutoff = 15/(Fs/2);
filtertype = 'low';
[b,a] = butter(filterorder,filtercutoff,filtertype);
clear Fs filterorder filtercutoff filtertype
% Right shin
Gyro_St_RS_r = filtfilt(b,a,Gyro_St_RS_r);
Gyro_St_LS_r = filtfilt(b,a,Gyro_St_LS_r);
Gyro_St_RS_x_r = filtfilt(b,a,Gyro_St_RS_x_r);
Acc_St_RS_z_r = filtfilt(b,a,Acc_St_RS_z_r);
Acc_St_RS_x_r = filtfilt(b,a,Acc_St_RS_x_r);
Acc_St_RS_y_r = filtfilt(b,a,Acc_St_RS_y_r);
[C,L] = wavedec(Acc_St_RS_x_r,4,'db5');
DG_x = wrcoef('a',C,L,'db5',4);
[C,L] = wavedec(Acc_St_RS_y_r,4,'db5');
DG_y = wrcoef('a',C,L,'db5',4);
accelometer = sqrt(DG_x .^ 2 + DG_y .^2)/2;
% Left shin
Gyro_St_LS_x_r = filtfilt(b,a,Gyro_St_LS_x_r);
Gyro_St_RS_y_r = filtfilt(b,a,Gyro_St_RS_y_r);
Gyro_St_LS_y_r = filtfilt(b,a,Gyro_St_LS_y_r);
Acc_St_LS_z_r = filtfilt(b,a,Acc_St_LS_z_r);
Acc_St_LS_x_r = filtfilt(b,a,Acc_St_LS_x_r);
Acc_St_LS_y_r = filtfilt(b,a,Acc_St_LS_y_r);
[C,L] = wavedec(Acc_St_LS_x_r,4,'db5');
DG_x_L = wrcoef('a',C,L,'db5',4);
[C,L] = wavedec(Acc_St_LS_y_r,4,'db5');
DG_y_L = wrcoef('a',C,L,'db5',4);
accelometer_L = sqrt(DG_x_L .^ 2 + DG_y_L .^2)/2;
clear a b
 
%% Define the propulsion phase by looking at the accelometer
Propulsion_point_R = FindPropulsion(accelometer,[res_St.HsR],[res_St.ToR]);
Propulsion_point_L = FindPropulsion(accelometer_L,[res_St(1:end-1).HsL],[res_St(2:end).ToL]);


%% plot the walking signal and define the Toe off and heel strike



pressure_Dir = 'Z:\Projects BCM\H-38994 BLANKET BCM\Studies\Fscan Frailty\Data\Raw Sensor data\Fscan Healthy';
cd(pressure_Dir);
[file, path,index] = uigetfile({'*.csv'},...
    'File Selector');
cd(path)
[pressure_Whole,~,~] =xlsread(file);

pressure = pressure_Whole(:,5);
figure('Position',[10,10,1800,1800]);
plot(pressure);

[points ~] = getpts;
points = points *2;
close;

% Test plot
[~, c1] = size(res_St);
for i = 1:c1-1
    if (file(8) == 'R' || file(7) == 'R')
        fig = figure(1);
        %                    h(1)=subplot(2,1,1); hold on;
        plot(res_St(i).HsR,Gyro_St_RS_r(res_St(i).HsR),'.k','markersize',12);
        %                    h(1)=subplot(2,1,1); hold on;
        plot(res_St(i).ToR,Gyro_St_RS_r(res_St(i).ToR),'*r','markersize',12);
        ylim([-300 600]);
    end
end

if (points > res_St(3).HsR)
    
    %                             h(1)=subplot(2,1,1); hold on;
    plot(pressure_resample_whole(round(points)-res_St(3).HsR:end) * 4,'g');
    % h(1)=subplot(2,1,1); hold on; plot(pressure_resample_Heel(round(points)-res_St(3).HsR:end),'black');
    [C,L] = wavedec(diff(pressure_resample_whole),5,'db5');
    A5 = wrcoef('a',C,L,'db5',4);
    %                             h(1)=subplot(2,1,1); hold on;
    plot(A5(round(points)-res_St(3).HsR:end) * 10,'black');
    
    [C,L] = wavedec(diff(diff(Gyro_St_RS_r)),4,'db5');
    DG = wrcoef('a',C,L,'db5',4);
    [C,L] = wavedec(diff(pressure_resample_whole),5,'db5');
    A5 = wrcoef('a',C,L,'db5',4);
    
    
    %  h(1)=subplot(2,1,1); hold on; plot(IntegGyro_St_RS_r,'k');
    %  h(1)=subplot(2,1,1); hold on; plot(DG* 50,'k');
    %  h(1)=subplot(1,1,1); hold on; plot(dpressureHeeldx(round(points)-res_St(3).HsR:end),'c');
    %                             h(1)=subplot(2,1,1); hold on;
    plot(A5(round(points)-res_St(3).HsR:end)* 10,'c');
    PressureNew = A5(round(points)-res_St(3).HsR:end);
    for u = 2:length(res_St)-1
        pointpressureMax = find(PressureNew(res_St(u).HsR + 40 :res_St(u).ToR) == max(PressureNew(res_St(u).HsR + 40:res_St(u).ToR)));
        %                                 h(1)=subplot(2,1,1); hold on;
        plot(res_St(u).HsR + pointpressureMax + 40,Gyro_St_RS_r(res_St(u).HsR + pointpressureMax + 40),'ok','markersize',12);
        
    end
else
    %                             h(1)=subplot(2,1,1); hold on;
    plot(res_St(3).HsR - round(points):length(pressure_resample_whole)+res_St(3).HsR - round(points)-1,pressure_resample_whole * 4,'g');
    % [C,L] = wavedec(diff(pressure_resample_whole),5,'db5');
    % A5 = wrcoef('a',C,L,'db5',4);
    % A5 = detrend(A5);
    % h(1)=subplot(2,1,1);  hold on; plot(res_St(3).HsR - round(points):length(pressure_resample_whole)+res_St(3).HsR - round(points)-1,pressure_resample_whole,'black');
    % h(1)=subplot(2,1,1); hold on; plot(IntegGyro_St_RS_r,'k');
    % h(1)=subplot(2,1,1); hold on; plot(Gyro_St_RS_x_r,'k');
    % h(1)=subplot(2,1,1); hold on;
    
    % [C,L] = wavedec(dpressureHeeldx,5,'db5');
    % a5 = wrcoef('a',C,L,'db5',4);
    [C,L] = wavedec(diff(diff(Gyro_St_RS_r)),4,'db5');
    DG = wrcoef('a',C,L,'db5',4);
    %                         h(1)=subplot(2,1,1); hold on; plot(DG * 50,'k');
    
    [C,L] = wavedec(diff(pressure_resample_whole),5,'db5');
    A5 = wrcoef('a',C,L,'db5',4);
    % plot(Gyro_St_RS_r)
    % plot(res_St(3).HsR - round(points):length(dpressureHeeldx)+res_St(3).HsR - round(points)-1,a5,'black');
    % h(1)=subplot(2,1,1); hold on;plot(res_St(3).HsR - round(points):length(A5)+res_St(3).HsR - round(points)-1,A5*10,'c');
    
    PressureNew = [zeros(res_St(3).HsR - round(points)-1,1); A5];
    
    
    for u = 2:length(res_St)-1
        pointpressureMax = find(PressureNew(res_St(u).HsR + 40:res_St(u).ToR) == max(PressureNew(res_St(u).HsR + 40:res_St(u).ToR)))
        %                                 h(1)=subplot(2,1,1); hold on;
        plot(res_St(u).HsR+pointpressureMax + 40,Gyro_St_RS_r(res_St(u).HsR +40 +pointpressureMax),'ok','markersize',12);
        
    end
end

%% Adding the propulsion to the Gait analysis
result = res_St;
for i = 1:size(res_St,2)
    result(i).Pro_R = Propulsion_point_R(1,i);
    if i ==1
        result(1).Pro_L = [];
    else
        if(Propulsion_point_L(1,i-1) ~= 1)
            result(i).Pro_L = Propulsion_point_L(1,i-1);
        else 
             result(i).Pro_L = [];
        end
    end
end

SavingDir = "Z:\Projects BCM\H-38994 BLANKET BCM\Studies\Fscan Frailty\Results\PropulsionPoints";

cd(SavingDir);
save(



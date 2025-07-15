clear all
close all
clc


RawDataDir  = 'Z:\Projects BCM\H-38994 BLANKET BCM\Studies\Fscan Frailty\Data\Raw Sensor data\Gait\';


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
        
%% 
dir_current = fullfile([directory listDir(s).name]);
[SensorData, samp_rate,~] = LoadLEGSysRawData5Sensors(dir_current);

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

Gyro_St_RS_r   = resample(Gyro_St_RS,200,samp_rate);
Gyro_St_RS_x_r = resample(Gyro_St_RS_x,200,samp_rate);
Gyro_St_RS_y_r = resample(Gyro_St_RS_y,200,samp_rate);
Acc_St_RS_z_r   = resample(Acc_St_RS_z,200,samp_rate);
Acc_St_RS_x_r = resample(Acc_St_RS_x,200,samp_rate);
Acc_St_RS_y_r = resample(Acc_St_RS_y,200,samp_rate);

Gyro_St_LS_r   = resample(Gyro_St_LS,200,samp_rate);
Gyro_St_LS_x_r = resample(Gyro_St_LS_x,200,samp_rate);
Gyro_St_LS_y_r = resample(Gyro_St_LS_y,200,samp_rate);
clear Gyro_St_RS Gyro_St_RS_x Gyro_St_RS_y
clear Gyro_St_LS Gyro_St_LS_x Gyro_St_LS_y

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
clear a b

% Test plot
fig = figure(1);
clf(fig)
h(1)=subplot(2,1,1); hold on; plot(Gyro_St_RS_r);
h(1)=subplot(2,1,1); hold on; plot(accelometer * 400);
%                     h(1)=subplot(2,1,1); hold on; plot(-DG_x * 400,'c');
%                      h(1)=subplot(2,1,1); hold on; plot(-DG_y * 400,'k');
xlim([1 length(Gyro_St_RS_r)])
grid on;

[start_points ~] = getpts
start_points = round(start_points);
%     close;

[end_points ~] = getpts
end_points = round(end_points);
max = 0;
s = 1;
num = 0;
num_positive = 0;
for i = start_points : 3 :end_points
    first_sample = [i:i+2];
    second_sample = [i+3:i+5];
    temp_mean1 = mean(accelometer(first_sample));
    temp_STD1 = std(accelometer(first_sample));
    temp_mean2 = mean(accelometer(second_sample));
    temp_STD2 = std(accelometer(second_sample));
    result(round(i/3)+1,:) = [i temp_mean2-temp_mean1/s temp_mean2-temp_mean1 temp_STD2-temp_STD1 temp_mean1 temp_mean2 temp_STD1 temp_STD2];
    if temp_mean2-temp_mean1 ~= 0
    s = temp_mean2-temp_mean1;
    end
    if s > .01
        num_positive = num_positive + 1;
    else num_positive =0;
    end
    if num_positive > 4
        point_interest = i - 9;
        m = FindMin(accelometer,point_interest)
        return;
    end
    
end


function m = FindMin(x,poin)
    [value, m] = min(x(poin-5:poin+5));
    m = poin - 5 + m;
end

                    
                    











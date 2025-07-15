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
RawDataDir  = 'Z:\Projects BCM\H-38994 BLANKET BCM\Studies\Fscan Frailty\Data\Converted\Gait';

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

PatientIndexData = s(1)
PatientName = listDir(PatientIndexData).name

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
clear PatientIndexData
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


%% Adding the propulsion to the Gait analysis
result = res_St;
for i = 1:size(res_St,2)
    result(i).Pro_R = Propulsion_point_R(1,i);
    if i ==1
        result(1).Pro_L  = [];
    else result(i).Pro_L = Propulsion_point_R(1,i);
    end
end

%% Analysis of the gyro paramters

%     cd('Z:\Personnel\Hadi\Projects\Fraility and Propulsion\Codes');
    
    Gyro_St_ss = [Gyro_St_RS_r, Gyro_St_LS_r];  
    res_St = GaitAnalyze(Gyro_St_ss, Height);
    [~, r] = size(res_St);
    res_St = res_St(2:r-1);
    
    angle_St_RS_x_r=(-integrationline(Gyro_St_RS_x_r)); % angle data for rotation around x axis
    angle_St_LS_x_r=(-integrationline(Gyro_St_LS_x_r));
    angle_St_RS_y_r=(-integrationline(Gyro_St_RS_y_r)); % angle data for rotation around y axis
    angle_St_LS_y_r=(-integrationline(Gyro_St_RS_y_r));
    
    abs_Gyro_RS_xy_r= sqrt(Gyro_St_RS_x_r.^2+Gyro_St_RS_y_r.^2); % Absoulte value for sum of angular velocities in xy plane
    abs_Gyro_LS_xy_r= sqrt(Gyro_St_LS_x_r.^2+Gyro_St_LS_y_r.^2);
    

%% calculating the parameters

% Redetection
    for n = 1:length(res_St)
        res_St(n).HsL;
        St_L_Hs(n,:)   = [Gyro_St_LS_r(res_St(n).HsL) res_St(n).HsL];
        St_L_To(n,:)   = [Gyro_St_LS_r(res_St(n).ToL) res_St(n).ToL];
        St_R_Hs(n,:)   = [Gyro_St_RS_r(res_St(n).HsR) res_St(n).HsR];
        St_R_To(n,:)   = [Gyro_St_RS_r(res_St(n).ToR) res_St(n).ToR];
        St_R_PSS(n)   = [res_St(n).PeakSwingSpeedR];                            % PSS: Peak Swing Speed
        St_L_PSS(n)   = [res_St(n).PeakSwingSpeedL];    
    end
    
    St_R_Hs(1,:)   = [];
    St_R_To(1,:)   = [];
    St_L_Hs(end,:) = [];
    St_L_To(1,:)   = [];
    St_R_PSS   =  St_R_PSS(2:end-1)';
    St_L_PSS   =  St_L_PSS(2:end-1)';

    %% Detection of stance phase
    St_StanceR = zeros(length(St_R_Hs),2);
    St_StanceL = zeros(length(St_L_Hs),2);
    
    peak_StanceR_x = zeros(length(St_R_Hs),2);                  % Peak RMS stance velocity in x-direction Gyroscopic data and location
    peak_StanceL_x = zeros(length(St_L_Hs),2);
    
    peak_SwingR_x = zeros(length(St_R_Hs)-1,2);
    peak_SwingL_x = zeros(length(St_L_Hs)-1,2);
    
    prop_max_acc_R = zeros(length(St_R_Hs),4);                   % Max propulsion acceleration and location
    prop_max_acc_L = zeros(length(St_L_Hs),4);
    
    peak_abs_xy_R = zeros(length(St_R_Hs),2);                   % Peaks of absoulte value for sum of angular velocities in xy plane
    peak_abs_xy_L = zeros(length(St_L_Hs),2);    
    peak_angle_R_x = zeros(length(St_R_Hs),2);                  % Peaks of angles of rotation around x axis  
    peak_angle_L_x = zeros(length(St_L_Hs),2);
    
    tro_angle_R_x = zeros(length(St_R_Hs),2);                  % Troughs of angles of rotation around x axis  
    tro_angle_L_x = zeros(length(St_L_Hs),2);
    
    peak_angle_R_y = zeros(length(St_R_Hs),2);                  % Peaks of angles of rotation around y axis  
    peak_angle_L_y = zeros(length(St_L_Hs),2);
    
    tro_angle_R_y = zeros(length(St_R_Hs),2);                  % Troughs of angles of rotation around y axis  
    tro_angle_L_y = zeros(length(St_L_Hs),2);
    
    for n = 1:length(St_R_Hs)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Right Side  
            [Rpks,Ridx] = findpeaks(Gyro_St_RS_r(St_R_Hs(n,2)+41:St_R_To(n,2)));
            [R_pk,R_id] = max(Rpks);
            St_StanceR(n,:) = [R_pk St_R_Hs(n,2)+40+Ridx(R_id)];

            prop_data_R =[(0.05*(St_R_Hs(n,2)+40+Ridx(R_id):St_R_To(n,2)))', Gyro_St_RS_r(St_R_Hs(n,2)+40+Ridx(R_id):St_R_To(n,2))];     % Propultion time and angular velocity values
            [~,~,R_max_acc,R_id_acc]= createFit(prop_data_R(:,1),prop_data_R(:,2));

            prop_max_acc_R(n,:)=[abs(R_max_acc),Gyro_St_RS_r(St_StanceR(n,2)+R_id_acc-1), St_StanceR(n,2)+R_id_acc-1,R_id_acc/length(prop_data_R(1,:))]; %updating the max propulsion acceleration (MPA), angular velocity value at MPA, MPA location and ratio of MPA/stance length

            % Finding max propulsion acceleration and location)
            [Rpks,Ridx] = findpeaks(rms(Gyro_St_RS_x_r(St_R_Hs(n,2):St_R_To(n,2)),2));  
            [R_pk,R_id] = max(Rpks);
            peak_StanceR_x(n,:) = [R_pk St_R_Hs(n,2)-1+Ridx(R_id)];

            m=n+1; 
            if m<length(St_R_Hs) % finding max angular velocity of swing in x direction
                [Rpks,Ridx] = findpeaks(rms(Gyro_St_RS_x_r(St_R_To(n,2):St_R_Hs(m,2)),2));
                [R_pk,R_id] = max(Rpks);
                peak_SwingR_x(n,:) = [R_pk St_R_To(n,2)-1+Ridx(R_id)];
            end

            [R_pk,R_id] = max(abs_Gyro_RS_xy_r(St_R_Hs(n,2)-200:St_R_Hs(n,2)+200)); % Finding the peaks of absoulte value for sum of angular velocities in xy plane
            peak_abs_xy_R (n,:) = [R_pk St_R_Hs(n,2)-201+R_id];

            [R_pk,R_id] = max(angle_St_RS_x_r(St_R_Hs(n,2)-200:St_R_To(n,2)+50)); % Finding peaks of angles of rotation around x axis
            peak_angle_R_x (n,:) = [R_pk St_R_Hs(n,2)-201+R_id];

            [R_pk,R_id] = min(angle_St_RS_x_r(St_R_Hs(n,2)-200:St_R_To(n,2)+50)); % Finding troughs of angles of rotation around x axis
            tro_angle_R_x (n,:) = [R_pk St_R_Hs(n,2)-201+R_id];

            [R_pk,R_id] = max(angle_St_RS_y_r(St_R_Hs(n,2)-200:St_R_To(n,2)+50)); % Finding peaks of angles of rotation around y axis
            peak_angle_R_y (n,:) = [R_pk St_R_Hs(n,2)-201+R_id];

            [R_pk,R_id] = min(angle_St_RS_y_r(St_R_Hs(n,2)-200:St_R_To(n,2)+50)); % Finding peaks of angles of rotation around y axis
            tro_angle_R_y (n,:) = [R_pk St_R_Hs(n,2)-201+R_id];        

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Left Side 
            [Lpks,Lidx] = findpeaks(Gyro_St_LS_r(St_L_Hs(n,2)+41:St_L_To(n,2)));
            [L_pk,L_id] = max(Lpks);
            St_StanceL(n,:) = [L_pk St_L_Hs(n,2)+40+Lidx(L_id)]
            
            L_max_acc = [];
            prop_data_L =[(0.05*(St_L_Hs(n,2)+40+Lidx(L_id):St_L_To(n,2)))', Gyro_St_LS_r(St_L_Hs(n,2)+40+Lidx(L_id):St_L_To(n,2))];
            [~,~,L_max_acc,L_id_acc]= createFit(prop_data_L(:,1),prop_data_L(:,2));
            
            
            if isempty(L_max_acc)
                continue;
            end
                
            prop_max_acc_L(n,:)=[abs(L_max_acc),Gyro_St_LS_r(St_StanceL(n,2)+L_id_acc-1), St_StanceL(n,2)+L_id_acc-1,L_id_acc/length(prop_data_L(1,:))];

            [Lpks,Lidx] = findpeaks(rms(Gyro_St_LS_x_r(St_L_Hs(n,2):St_L_To(n,2)),2));
            [L_pk,L_id] = max(Lpks);
            peak_StanceL_x(n,:) = [L_pk St_L_Hs(n,2)-1+Lidx(L_id)];

            if m < length(St_R_Hs)
                [Lpks,Lidx] = findpeaks(rms(Gyro_St_LS_x_r(St_L_To(n,2):St_L_Hs(m,2)),2));
                [L_pk,L_id] = max(Lpks);
                peak_SwingL_x(n,:) = [L_pk St_L_To(n,2)-1+Lidx(L_id)];
            end

            [L_pk,L_id] = max(abs_Gyro_LS_xy_r(St_L_Hs(n,2)-200:St_L_Hs(n,2)+200));
            peak_abs_xy_L (n,:) = [L_pk St_L_Hs(n,2)-201+L_id];

            [L_pk,L_id] = max(angle_St_LS_x_r(St_L_Hs(n,2)-200:St_L_To(n,2)+50)); 
            peak_angle_L_x (n,:) = [L_pk St_L_Hs(n,2)-201+L_id];

            [L_pk,L_id] = min(angle_St_LS_x_r(St_L_Hs(n,2)-200:St_L_To(n,2)+50)); 
            tro_angle_L_x (n,:) = [L_pk St_L_Hs(n,2)-201+L_id];

            [L_pk,L_id] = max(angle_St_LS_y_r(St_L_Hs(n,2)-200:St_L_To(n,2)+50)); 
            peak_angle_L_y (n,:) = [L_pk St_L_Hs(n,2)-201+L_id];

            [L_pk,L_id] = min(angle_St_LS_y_r(St_L_Hs(n,2)-200:St_L_To(n,2)+50)); 
            tro_angle_L_y (n,:) = [L_pk St_L_Hs(n,2)-201+L_id];
    end
    %% Calculation of parameters
    Param_St_L_Hs_Pk = St_StanceL - St_L_Hs;
    Param_St_L_Pk_To = St_L_To    - St_StanceL;
    Param_St_L_Hs_To = St_L_To    - St_L_Hs;
    Param_St_R_Hs_Pk = St_StanceR - St_R_Hs;
    Param_St_R_Pk_To = St_R_To    - St_StanceR;
    Param_St_R_Hs_To = St_R_To    - St_R_Hs;
    
    xParam_ST_R_Pk_Pk = prop_max_acc_R(:,[2 3])- peak_StanceR_x;
    xParam_ST_L_Pk_Pk = prop_max_acc_L(:,[2 3]) - peak_StanceL_x;
    
    Param_x_R_Pk_Tr = abs(peak_angle_R_x-tro_angle_R_x);
    Param_y_R_Pk_Tr = abs(peak_angle_R_y-tro_angle_R_y);
    
    Param_x_L_Pk_Tr = abs(peak_angle_L_x-tro_angle_L_x);
    Param_y_L_Pk_Tr = abs(peak_angle_L_y-tro_angle_L_y);
    
    Param_St_L_Hs_Pk(:,2) = Param_St_L_Hs_Pk(:,2)/2*10; % unit: ms
    Param_St_L_Pk_To(:,2) = Param_St_L_Pk_To(:,2)/2*10;
    Param_St_L_Hs_To(:,2) = Param_St_L_Hs_To(:,2)/2*10;
    Param_St_R_Hs_Pk(:,2) = Param_St_R_Hs_Pk(:,2)/2*10;
    Param_St_R_Pk_To(:,2) = Param_St_R_Pk_To(:,2)/2*10;
    Param_St_R_Hs_To(:,2) = Param_St_R_Hs_To(:,2)/2*10;
    
    xParam_ST_R_Pk_Pk(:,2) = xParam_ST_R_Pk_Pk(:,2)/2*10;
    xParam_ST_L_Pk_Pk(:,2) = xParam_ST_L_Pk_Pk(:,2)/2*10;
    
    Param_x_R_Pk_Tr(:,2) = Param_x_R_Pk_Tr(:,2)/2*10;
    Param_y_R_Pk_Tr(:,2) = Param_y_R_Pk_Tr(:,2)/2*10;
    Param_x_L_Pk_Tr(:,2) = Param_x_L_Pk_Tr(:,2)/2*10;
    Param_y_L_Pk_Tr(:,2) = Param_y_L_Pk_Tr(:,2)/2*10;
    
    if St_L_To(1,2) < St_R_To(1,2)
        Param_St_L_Hs_R_Hs = St_R_Hs          - St_L_Hs;
        Param_St_R_Hs_L_Hs = St_L_Hs(2:end,:) - St_R_Hs(1:end-1,:);
    else
        Param_St_L_Hs_R_Hs = St_R_Hs(2:end,:) - St_L_Hs(1:end-1,:);
        Param_St_R_Hs_L_Hs = St_L_Hs          - St_R_Hs;
    end
    Result_St_L = [Param_St_L_Hs_Pk Param_St_L_Pk_To Param_St_L_Hs_To xParam_ST_L_Pk_Pk Param_x_L_Pk_Tr Param_y_L_Pk_Tr];
    Result_St_R = [Param_St_R_Hs_Pk Param_St_R_Pk_To Param_St_R_Hs_To xParam_ST_R_Pk_Pk Param_x_R_Pk_Tr Param_y_R_Pk_Tr];
    
    slope_R = -Result_St_R(:,3)./Result_St_R(:,4)*1000;     % Average slope of angular velocity for midstance peak to toe-off
    slope_L = -Result_St_L(:,3)./Result_St_L(:,4)*1000;
    prop_time_ratio_R = Result_St_R(:,4)./Result_St_R(:,6); % Propultion to stance time ratio
    prop_time_ratio_L = Result_St_L(:,4)./Result_St_L(:,6);
    
    slope_R_angle_x = Result_St_R(:,9)./Result_St_R(:,10)*1000;     % Average slope of angle x-axis
    slope_L_angle_x = Result_St_L(:,9)./Result_St_L(:,10)*1000;
    
    slope_R_angle_y = Result_St_R(:,11)./Result_St_R(:,12)*1000;     % Average slope of angle y-axis
    slope_L_angle_y = Result_St_L(:,11)./Result_St_L(:,12)*1000;
    
    
    % Output
        Outresults(1).ID = PatientName;
        % Peak Swing Speed outputs (PSS)
            Outresults(1).Peak_swing_speed_outputs_R_Mean = mean(St_R_PSS);
            Outresults(1).Peak_swing_speed_outputs_R_Std  = std(St_R_PSS);
            Outresults(1).Peak_swing_speed_outputs_R_Cov  = std(St_R_PSS)/mean(St_R_PSS);
            Outresults(1).Peak_swing_speed_outputs_L_Mean = mean(St_L_PSS);
            Outresults(1).Peak_swing_speed_outputs_L_Std  = std(St_L_PSS);
            Outresults(1).Peak_swing_speed_outputs_L_Cov  = std(St_L_PSS)/mean(St_L_PSS);
         % Max Propulsion acceleration
            Outresults(1).Max_Propulsion_acceleration_R_Mean = mean(prop_max_acc_R(:,1));
            Outresults(1).Max_Propulsion_acceleration_R_Std  = std(prop_max_acc_R(:,1));
            Outresults(1).Max_Propulsion_acceleration_R_Cov  = std(prop_max_acc_R(:,1))/mean(prop_max_acc_R(:,1));
            Outresults(1).Max_Propulsion_acceleration_L_Mean = mean(prop_max_acc_L(:,1));
            Outresults(1).Max_Propulsion_acceleration_L_Std  = std(prop_max_acc_L(:,1));
            Outresults(1).Max_Propulsion_acceleration_L_Cov  = std(prop_max_acc_L(:,1))/mean(prop_max_acc_L(:,1));
%%
    SavePath = 'Z:\Projects BCM\H-38994 BLANKET BCM\Studies\Fscan Frailty\Results\Patients\'
    cd(SavePath);
    str = date;
    Output = [PatientName,'_Gait_Results_',str,'.xlsx'];
    if  exist(Output,'file')==2
        delete(Output)
    end
    writetable(struct2table(Outresults), [SavePath,Output])
    cd(CodeDir)
    display(['Results Saved (Excel): ', Output]);
    clear ResultOutput  CurDir  DataPath
    return;
    
    
    
   output_results(zz,:)=[sep_, ...
        mean(-Result_St_R(:,3)),std(-Result_St_R(:,3)),std(-Result_St_R(:,3))/mean(-Result_St_R(:,3)),...           % Midstance peak velocity outputs (delta V) (PMS)
        mean(-Result_St_L(:,3)),std(-Result_St_L(:,3)),std(-Result_St_L(:,3))/mean(-Result_St_L(:,3)),...
        mean(slope_R),std(slope_R),std(slope_R)/mean(slope_R),...                                                    % Average slope of angular velocity for midstance peak to toe-off outputs
        mean(slope_L),std(slope_L),std(slope_L)/mean(slope_L),...
        mean(prop_time_ratio_R),std(prop_time_ratio_R),std(prop_time_ratio_R)/mean(prop_time_ratio_R),...            % Propultion/stance time ratio
        mean(prop_time_ratio_L),std(prop_time_ratio_L),std(prop_time_ratio_L)/mean(prop_time_ratio_L),...
        mean(Result_St_R(:,6)),std(Result_St_R(:,6)),std(Result_St_R(:,6))/mean(Result_St_R(:,6)), ...                   % Absolute Propulsion time (delta T)
        mean(Result_St_L(:,6)),std(Result_St_L(:,6)),std(Result_St_L(:,6))/mean(Result_St_L(:,6)),...
        mean(prop_max_acc_R(:,1)),std(prop_max_acc_R(:,1)),std(prop_max_acc_R(:,1))/mean(prop_max_acc_R(:,1)), ...       % Max Propulsion acceleration
        mean(prop_max_acc_L(:,1)),std(prop_max_acc_L(:,1)),std(prop_max_acc_L(:,1))/mean(prop_max_acc_L(:,1)),...
        mean(prop_max_acc_R(:,4)),std(prop_max_acc_R(:,4)),std(prop_max_acc_R(:,4))/mean(prop_max_acc_R(:,4)),...        % Percentile of location of max Propulsion Acceleration within the propulsion
        mean(prop_max_acc_L(:,4)),std(prop_max_acc_L(:,4)),std(prop_max_acc_L(:,4))/mean(prop_max_acc_L(:,4)),...
        mean(-Result_St_R(2:end,3)./St_R_PSS),std(-Result_St_R(2:end,3)./St_R_PSS),std(-Result_St_R(2:end,3)./St_R_PSS)/mean(-Result_St_R(2:end,3)./St_R_PSS),...      % Ratio of PMS/PSS
        mean(-Result_St_L(2:end,3)./St_L_PSS),std(-Result_St_L(2:end,3)./St_L_PSS),std(-Result_St_L(2:end,3)./St_L_PSS)/mean(-Result_St_L(2:end,3)./St_L_PSS),...
        FFT_Max_R,FFT_Max_L,...                                                                                           % max frequecy
        FFT_Max_R_x,FFT_Max_R_x,...
        mean(peak_StanceR_x(:,1)),std(peak_StanceR_x(:,1)),std(peak_StanceR_x(:,1))/mean(peak_StanceR_x(:,1)), ...          % Peak stance velocity in x-direction Gyroscopic data
        mean(peak_StanceL_x(:,1)),std(peak_StanceL_x(:,1)),std(peak_StanceL_x(:,1))/mean(peak_StanceL_x(:,1)), ...
        mean(Result_St_R(:,8)),std(Result_St_R(:,8)),std(Result_St_R(:,8))/mean(Result_St_R(:,8)), ...                   % Time difference between the peak swing (z) and peak stance (x)
        mean(Result_St_L(:,8)),std(Result_St_L(:,8)),std(Result_St_L(:,8))/mean(Result_St_L(:,8)),...
        mean(peak_SwingR_x(:,1)),std(peak_SwingR_x(:,1)),std(peak_SwingR_x(:,1))/mean(peak_SwingR_x(:,1)), ...          % Peak swing velocity in x-axis Gyroscopic data
        mean(peak_SwingL_x(:,1)),std(peak_SwingL_x(:,1)),std(peak_SwingL_x(:,1))/mean(peak_SwingL_x(:,1)), ...
        mean(Result_St_R(:,9)),std(Result_St_R(:,9)),std(Result_St_R(:,9))/mean(Result_St_R(:,9)), ...                   % Change in angle per gait cycle x-axis
        mean(Result_St_L(:,9)),std(Result_St_L(:,9)),std(Result_St_L(:,9))/mean(Result_St_L(:,9)),...
        mean(Result_St_R(:,11)),std(Result_St_R(:,11)),std(Result_St_R(:,11))/mean(Result_St_R(:,11)), ...                   % Change in angle per gait cycle y-axis
        mean(Result_St_L(:,11)),std(Result_St_L(:,11)),std(Result_St_L(:,11))/mean(Result_St_L(:,11)),...
        mean(slope_R_angle_x),std(slope_R_angle_x),std(slope_R_angle_x)/mean(slope_R_angle_x),...                   % Slope of angle per gait cycle x-axis
        mean(slope_L_angle_x),std(slope_L_angle_x),std(slope_L_angle_x)/mean(slope_L_angle_x),...
        mean(slope_R_angle_y),std(slope_R_angle_y),std(slope_R_angle_y)/mean(slope_R_angle_y),...                    % Slope of angle per gait cycle y-axis
        mean(slope_L_angle_y),std(slope_L_angle_y),std(slope_L_angle_y)/mean(slope_L_angle_y),...
        mean(peak_abs_xy_R(:,1)),std(peak_abs_xy_R(:,1)),std(peak_abs_xy_R(:,1))/mean(peak_abs_xy_R(:,1)),...         % Peak absolute xy angular velocity
        mean(peak_abs_xy_L(:,1)),std(peak_abs_xy_L(:,1)),std(peak_abs_xy_L(:,1))/mean(peak_abs_xy_L(:,1)),...
        ];
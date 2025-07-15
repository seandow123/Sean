%--------------------------------------------------------------------------
% Project Frailty and Propulsion
%   The purpose of the function is to calculate the propulsion data using
%   the legsys sensor
% Orginal Code Design : Hyoki Lee (Last Version 1_2, see Previous Version)
% Author:  Hung Nguyen
% Code Description:
%   Version 1.0
%   This code will load in the legsys sensor dat and use the algorithm to
%   detect the propulsion parameters. The propulsion parameter will follow
%   the same naming scheme as the propulsion paper.
%--------------------------------------------------------------------------
%% Initialization
clear all; clc; home;
warning off
option_sensor   = 1;          % 1: Use only Wrist Sensor 2: Use Wrist + Upper Arm Sensor
format short g
%% Options for saving parameter results,number of sensor and executed task.
%     option_save = 1;          % 1: Save the parameter results as Excel sheet, otherwise: Not saved
%     option_task = 'STW'; 
    DataSet     = 'OHI ';      % Name of data set to analyzed ['CIPN' or 'PAD']

    % Get Version of the code
        CodeName = mfilename
        Versionindex    = max(find(CodeName=='v' | CodeName=='V'));
        CodeVersion     = CodeName(Versionindex:end);
        clear Versionindex
%% Load Starting point
%{
    Load in the excel files that contains the subject numbers and the
    begining and starting point.  The matrix is organized as
    [subjectNumber BegPoint EndPoint]
%}
    CurDir = pwd;
    ExcelDir = 'P:\Projects BCM\H-38994 BLANKET BCM\Studies\Propulsion\Results\';
    cd(ExcelDir)
    [PtsRange,~,~] =xlsread('OHi_Points.xlsx');
    
%% Directory Location
    % Code Directory
    CodeDir = fullfile(pwd);
    % Raw data directory
    RawDataDir  = 'P:\Projects BCM\H-38994 BLANKET BCM\Studies\Propulsion\Data\Raw Sensor data\Gait\';
    % Directory to same matlab and excel resutls
    ResultsDir  = 'P:\Projects BCM\H-38994 BLANKET BCM\Studies\Propulsion\Results\';
    
    %%  Load Data for Heel Area
    
    Area_Excel_Dir = 'P:\Projects BCM\H-38994 BLANKET BCM\Studies\Propulsion\Data\Analyzed\Fscan OHI';
    cd(Area_Excel_Dir);
    [Area_Heel,~,~] =xlsread('OHI 001L 6M _ heel_G.csv');
    [pressure_Heel,~,~] =xlsread('OHI 016R 6M_Heel_G_Contact pressure_Heel.csv');
    [pressure_Whole,~,~] =xlsread('OHI 016R 6M_Heel_G_Contact pressure_Whole.csv');
%     force = Area_Heel(:,5);
%     figure('Position',[10,10,1800,1800]);
%     plot(force);
%     [points ~] = getpts;
%     points = points *2;
%     close;
    
    pressure = pressure_Heel(:,5);
    figure('Position',[10,10,1800,1800]);
    plot(pressure);
    
%     return;
    
    [points ~] = getpts;
    points = points *2;
    close;
    
    
    
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
for zz = 1:length(s)        % length of selected folder for analysis
    % Selected Patient ID ( should be in the variable s:selection)
        % to access the direction use: listdir(index).name
        % patient_index based on the converted data directory(not patient
        % information
        PatientIndexData = s(zz);
        
        % patient name (OHI XXX)
        LenName    = length(DataSet);
        LenNameMax = LenName + 3;
        PatientNameData = listDir(PatientIndexData).name
        
    % Check to see if there are remark on the data - shorten the name
        len = length(PatientNameData);
        if len > LenNameMax  % If the folder name includes remarks                                                                                    
            PatientNameData = PatientNameData(1:find(PatientNameData==' ',1)-1);
        end
        
    % Converting to the real number of patient from string filename
        BegPts        = LenName + 1; % begin position of number. Ex. OHI 001 - begin of 001
        PatientNameDataDble = str2double(PatientNameData(BegPts:end));  
        PatientNameIndex = find(PtsRange(:,1)==PatientNameDataDble);
        clear BegPts len LenNameMax LenName 
        
    % Get start and end points of each subject
        StartPt = PtsRange(PatientNameIndex,2);
        EndPt   = PtsRange(PatientNameIndex,3);
        Height  = PtsRange(PatientNameIndex,4);
        
        if( isnan(StartPt) | isnan(EndPt))
            error('Invalid start/end point for subject')
        end
        
    % Load Sensor data
        %{
            * For propulsion Code (loadLegSysRawDataV3) there is quaternion
              output
            * SensorData = [s1 s2 s3 s4 s5]
                each sensor value is 10-tuple:
                1: time stamp
                [2 3 4]    : Gyrox Gyroy Gyroz
                [5 6 7]    : Accx  Accy  Accz
                [8 9 10 11]: w x y z
            * samp_rate = sample rate in Hz
            * pathname = current data directory
        %}
        dir_current = fullfile([directory listDir(PatientIndexData).name],filesep);
        clear PatientIndexData
        disp([' '])
        disp(['----- Analyzing ', PatientNameData,' ------'])
        
        % Intialize 
        SensorData = {};
        [SensorData, samp_rate,~] = LoadLEGSysRawData5Sensors(dir_current);
        if(isempty(SensorData))
            continue
        end
        if (strcmp(PatientNameData,'OHI 008'))
            SensorData = -SensorData;
        end
        
   
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
            
            if(SensorSetup == 2 || SensorSetup ==5)
                % Sensor 1 or "Right shin"
                Gyro_St_RS(:,1) =  SensorData(StartPt:EndPt,2);	% x for s1
                Gyro_St_RS(:,2) =  SensorData(StartPt:EndPt,3);	% y for s1
                Gyro_St_RS(:,3) = -SensorData(StartPt:EndPt,4);	% z for s1
                
                % Sensor 2 or "Left Shin"
                Gyro_St_LS(:,1) =  SensorData(StartPt:EndPt,22);	% x for s3
                Gyro_St_LS(:,2) =  SensorData(StartPt:EndPt,23);	% y for s3
                Gyro_St_LS(:,3) = -SensorData(StartPt:EndPt,24);	% z for s3
            else
                error('Wrong sensor setup')  
            end

        % Test plot
            x1 = 1;
            if (x1 == 1)
                x2 = max(length(Gyro_St_RS(:,3)));
            else
                x2 = 1000 + x1;
            end

            fig = figure(1);
            clf(fig);
            h(1)=subplot(2,1,1); plot(Gyro_St_RS(:,3));
            title([PatientNameData,' Right Shin Gyro about z'])
            grid on;
            xlim([x1 x2])

            h(2)=subplot(2,1,2); plot(Gyro_St_LS(:,3));
            title([PatientNameData,' Left Shin Gyro about z'])
            xlim([x1 x2])
            grid on
            linkaxes(h,'x')
            clear Gyro_St_RS Gyro_St_LS x1 x2 fig h  
            
        % Filter data
            % Right shin
                Gyro_St_RS_x =  filtfilt([1,-1],[1,-0.995],-SensorData(StartPt:EndPt,2));
                Gyro_St_RS_y =  filtfilt([1,-1],[1,-0.995],-SensorData(StartPt:EndPt,3));
                Gyro_St_RS   =  filtfilt([1,-1],[1,-0.995],-SensorData(StartPt:EndPt,4));	% z
            % Left shin
                Gyro_St_LS_x =  filtfilt([1,-1],[1,-0.995],-SensorData(StartPt:EndPt,22));
                Gyro_St_LS_y =  filtfilt([1,-1],[1,-0.995],-SensorData(StartPt:EndPt,23));
                Gyro_St_LS   =  filtfilt([1,-1],[1,-0.995],-SensorData(StartPt:EndPt,24));  % z
                clear StartPt EndPt
                
                % Test plot
                    fig = figure(1);
                    h(1)=subplot(2,1,1); hold on; plot(Gyro_St_RS);
                    ylabel('Angular Velocity')

                    h(2)=subplot(2,1,2); hold on; plot(Gyro_St_LS);
                    grid on
                    linkaxes(h,'x')
                    ylabel('Angular Velocity')
                    clear fig h  
        % Resampling
            Gyro_St_RS_r   = resample(Gyro_St_RS,200,samp_rate);
            Gyro_St_RS_x_r = resample(Gyro_St_RS_x,200,samp_rate);
            Gyro_St_RS_y_r = resample(Gyro_St_RS_y,200,samp_rate);

            Gyro_St_LS_r   = resample(Gyro_St_LS,200,samp_rate);
            Gyro_St_LS_x_r = resample(Gyro_St_LS_x,200,samp_rate);
            Gyro_St_LS_y_r = resample(Gyro_St_LS_y,200,samp_rate);
            clear Gyro_St_RS Gyro_St_RS_x Gyro_St_RS_y
            clear Gyro_St_LS Gyro_St_LS_x Gyro_St_LS_y
            
                % Test plot
                    fig = figure(1);
                    clf(fig)
                    h(1)=subplot(2,1,1); plot(Gyro_St_RS_r);
                    xlim([1 length(Gyro_St_RS_r)])
                    grid on;

                    h(2)=subplot(2,1,2);plot(Gyro_St_LS_r);
                    xlim([1 length(Gyro_St_RS_r)])
                    
                    grid on
                    linkaxes(h,'xy')
                    clear fig h  
        
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
            % Left shin    
                Gyro_St_LS_x_r = filtfilt(b,a,Gyro_St_LS_x_r);
                Gyro_St_RS_y_r = filtfilt(b,a,Gyro_St_RS_y_r);
                Gyro_St_LS_y_r = filtfilt(b,a,Gyro_St_LS_y_r);
                clear a b
                
                    % Test plot
                    fig = figure(1);
                    clf(fig)
                    h(1)=subplot(2,1,1); hold on; plot(Gyro_St_RS_r);
                    xlim([1 length(Gyro_St_RS_r)])
                    grid on;

                    h(2)=subplot(2,1,2);hold on; plot(Gyro_St_LS_r);
                    xlim([1 length(Gyro_St_RS_r)])
                    grid on;
                    
                    linkaxes(h,'xy')
                    clear fig h                  
    % Analysis
        Gyro_St_ss = [Gyro_St_RS_r, Gyro_St_LS_r];
        res_St     = GaitAnalyze(Gyro_St_ss, Height);
        clear Height 
        
        
        
        
           % Test plot
           [~, c1] = size(res_St)
           for i =1:c1
                    fig = figure(1);
                    h(1)=subplot(2,1,1); hold on; plot(res_St(i).HsR,Gyro_St_RS_r(res_St(i).HsR),'.k','markersize',12);
                    h(1)=subplot(2,1,1); hold on; plot(res_St(i).ToR,Gyro_St_RS_r(res_St(i).ToR),'*r','markersize',12);

%                     xlim([1 length(Gyro_St_RS_r)])
%                     grid on;

                    h(2)=subplot(2,1,2); hold on; plot(res_St(i).HsL,Gyro_St_LS_r(res_St(i).HsL),'.k','markersize',12);
                    h(2)=subplot(2,1,2); hold on; plot(res_St(i).ToL,Gyro_St_LS_r(res_St(i).ToL),'*r','markersize',12);
%                     area_resample_heel = resample(Area_Heel(:,5),2,1) * 10;
                    pressure_resample_heel = resample(pressure_Heel(:,5),2,1);
                    pressure_resample_whole = resample(pressure_Whole(:,5),2,1);
                    xlim([1 length(pressure_resample_heel)])
                    ylim([-300 700])
%                     grid on;
                    
                    linkaxes(h,'xy')
%                     clear fig h      
           end
%            if (points > res_St(2).HsR)
%                h(1)=subplot(2,1,1); hold on; plot(area_resample_heel(round(points)-res_St(2).HsR:end));  
%            else
%                h(1)=subplot(2,1,1); hold on; plot(res_St(2).HsR - round(points):length(area_resample_heel)+res_St(2).HsR - round(points)-1,area_resample_heel);
%            end

            if (points > res_St(2).HsR)
               h(1)=subplot(2,1,1); hold on; plot(pressure_resample_heel(round(points)-res_St(2).HsR:end));  
               h(1)=subplot(2,1,1); hold on; plot(pressure_resample_whole(round(points)-res_St(2).HsR:end)); 
           else
               h(1)=subplot(2,1,1); hold on; plot(res_St(2).HsR - round(points):length(pressure_resample_heel)+res_St(2).HsR - round(points)-1,pressure_resample_heel);
               h(1)=subplot(2,1,1); hold on; plot(res_St(2).HsR - round(points):length(pressure_resample_whole)+res_St(2).HsR - round(points)-1,pressure_resample_whole);
           end
                    
           return; 
                    
        
        
        
        % Angle data for rotation around x axis
        angle_St_RS_x_r =(-integrationline(Gyro_St_RS_x_r)); 
        angle_St_LS_x_r =(-integrationline(Gyro_St_LS_x_r));
        
        % angle data for rotation around y axis
        angle_St_RS_y_r =(-integrationline(Gyro_St_RS_y_r)); 
        angle_St_LS_y_r =(-integrationline(Gyro_St_LS_y_r));
        
        % Absoulte value for sum of angular velocities in xy plane
        abs_Gyro_RS_xy_r = sqrt(Gyro_St_RS_x_r.^2+Gyro_St_RS_y_r.^2); 
        abs_Gyro_LS_xy_r = sqrt(Gyro_St_LS_x_r.^2+Gyro_St_LS_y_r.^2);
                
        
         % Max frequecy extraction
            % Right Shin
                FFTn = 1024;
                FFT_R = abs(fft(isol_signal(Gyro_St_RS_r),FFTn));
                FFT_idx_R = (0:FFTn-1)/FFTn*200;
                [~,FFT_max_R] = max(medfilt1(FFT_R,10));
                FFT_Max_R = FFT_idx_R(FFT_max_R);
                clear FFT_idx_R
                

                FFT_R_x = abs(fft(isol_signal(Gyro_St_RS_x_r),FFTn));
                FFT_idx_R_x = (0:FFTn-1)/FFTn*200;
                [~,FFT_max_R_x] = max(medfilt1(FFT_R_x,10));
                FFT_Max_R_x = FFT_idx_R_x(FFT_max_R_x);
                
            % Left shin
                FFT_L = abs(fft(isol_signal(Gyro_St_LS_r),FFTn));
                FFT_idx_L = (0:FFTn-1)/FFTn*200;
                [~,FFT_max_L] = max(medfilt1(FFT_L,10));
                FFT_Max_L = FFT_idx_L(FFT_max_L);
                clear FFT_idx_R


                FFT_L_x = abs(fft(isol_signal(Gyro_St_LS_x_r),FFTn));
                FFT_idx_L_x = (0:FFTn-1)/FFTn*200;
                [~,FFT_max_L_x] = max(medfilt1(FFT_L_x,10));
                FFT_Max_L_x = FFT_idx_L_x(FFT_max_L_x);
                clear FFTn FFT_R FFT_L
            
        % Re-Detection
            [r1,c1] = size(res_St);
            clear r1 c1
            if (PatientNameDataDble==45)
                for n = 1:length(res_St)-4
                    St_L_Hs(n,:) = [Gyro_St_LS_r(res_St(n+3).HsL) res_St(n+3).HsL];
                    St_L_To(n,:) = [Gyro_St_LS_r(res_St(n+3).ToL) res_St(n+3).ToL];
                    St_R_Hs(n,:) = [Gyro_St_RS_r(res_St(n+3).HsR) res_St(n+3).HsR];
                    St_R_To(n,:) = [Gyro_St_RS_r(res_St(n+3).ToR) res_St(n+3).ToR];

                    % PSS: Peak Swing Speed
                    St_R_PSS(n)  = [res_St(n+3).PeakSwingSpeedR];                            
                    St_L_PSS(n)  = [res_St(n+3).PeakSwingSpeedL];
                end
            
            else    
                for n = 1:length(res_St)-2
                    St_L_Hs(n,:) = [Gyro_St_LS_r(res_St(n+1).HsL) res_St(n+1).HsL];
                    St_L_To(n,:) = [Gyro_St_LS_r(res_St(n+1).ToL) res_St(n+1).ToL];
                    St_R_Hs(n,:) = [Gyro_St_RS_r(res_St(n+1).HsR) res_St(n+1).HsR];
                    St_R_To(n,:) = [Gyro_St_RS_r(res_St(n+1).ToR) res_St(n+1).ToR];

                    % PSS: Peak Swing Speed
                    St_R_PSS(n)  = [res_St(n+1).PeakSwingSpeedR];                            
                    St_L_PSS(n)  = [res_St(n+1).PeakSwingSpeedL];
                end
            end
            
            
            St_R_Hs(1,:)   = [];
            St_R_To(1,:)   = [];
            St_L_Hs(end,:) = [];
            St_L_To(1,:)   = [];
            St_R_PSS   =  St_R_PSS(2:end-1)';
            St_L_PSS   =  St_L_PSS(2:end-1)';
        
                % Detection of stance phase
                St_StanceR = zeros(length(St_R_Hs),2);
                St_StanceL = zeros(length(St_L_Hs),2);

                % Peak RMS stance velocity in x-direction Gyroscopic data and location
                peak_StanceR_x = zeros(length(St_R_Hs),2);                  
                peak_StanceL_x = zeros(length(St_L_Hs),2);

                peak_SwingR_x = zeros(length(St_R_Hs)-1,2);
                peak_SwingL_x = zeros(length(St_L_Hs)-1,2);

                % Max propulsion acceleration and location
                prop_max_acc_R = zeros(length(St_R_Hs),4);                  
                prop_max_acc_L = zeros(length(St_L_Hs),4);

                % Peaks of absoulte value for sum of angular velocities in xy plane
                peak_abs_xy_R = zeros(length(St_R_Hs),2);                   
                peak_abs_xy_L = zeros(length(St_L_Hs),2);    

                % Peaks of angles of rotation around x axis  
                peak_angle_R_x = zeros(length(St_R_Hs),2);                  
                peak_angle_L_x = zeros(length(St_L_Hs),2);

                % Troughs of angles of rotation around x axis  
                tro_angle_R_x = zeros(length(St_R_Hs),2);                 
                tro_angle_L_x = zeros(length(St_L_Hs),2);

                % Peaks of angles of rotation around y axis  
                peak_angle_R_y = zeros(length(St_R_Hs),2);                  
                peak_angle_L_y = zeros(length(St_L_Hs),2);

                % Troughs of angles of rotation around y axis
                tro_angle_R_y = zeros(length(St_R_Hs),2);                    
                tro_angle_L_y = zeros(length(St_L_Hs),2);
                
                mm = length(St_R_Hs)
                P1=1;
                P2 = mm-1;
                if (PatientNameDataDble==13)
                    P2 = mm-2;
                elseif (PatientNameDataDble==14)
                    P2 = 7;
                elseif (PatientNameDataDble==16)
                    P1=2;
                elseif (PatientNameDataDble==43)
                    P1=2;
                elseif (PatientNameDataDble==44)
                    P2=4;
                elseif (PatientNameDataDble==451)
                    P1=3;
                end
                clear mm
                
                for n = P1:P2
                    n;
                    
                    % Right Side  
                        [Rpks,Ridx] = findpeaks(Gyro_St_RS_r(St_R_Hs(n,2)+41:St_R_To(n,2)));
                        [R_pk,R_id] = max(Rpks);
                        St_StanceR(n,:) = [R_pk St_R_Hs(n,2)+40+Ridx(R_id)];
                    
                        % Propultion time and angular velocity values
                        prop_data_R =[(0.05*(St_R_Hs(n,2)+40+Ridx(R_id):St_R_To(n,2)))', Gyro_St_RS_r(St_R_Hs(n,2)+40+Ridx(R_id):St_R_To(n,2))];     
                        [~,~,R_max_acc,R_id_acc]= createFit(prop_data_R(:,1),prop_data_R(:,2));
                        
                        
                        
                        %updating the max propulsion acceleration (MPA), angular velocity value at MPA, MPA location and ratio of MPA/stance length
                        prop_max_acc_R(n,:)=[abs(R_max_acc),Gyro_St_RS_r(St_StanceR(n,2)+R_id_acc-1), St_StanceR(n,2)+R_id_acc-1,R_id_acc/length(prop_data_R(1,:))]; 

                        % Finding max propulsion acceleration and location)
                        [Rpks,Ridx] = findpeaks(rms(Gyro_St_RS_x_r(St_R_Hs(n,2):St_R_To(n,2)),2));  
                        [R_pk,R_id] = max(Rpks);
                        peak_StanceR_x(n,:) = [R_pk St_R_Hs(n,2)-1+Ridx(R_id)]
        
                        m=n+1;
                        % finding max angular velocity of swing in x direction
                            if m<length(St_R_Hs) 
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
                        
                    % Left Side 
                        [Lpks,Lidx] = findpeaks(Gyro_St_LS_r(St_L_Hs(n,2)+41:St_L_To(n,2)));
                        [L_pk,L_id] = max(Lpks);
                        St_StanceL(n,:) = [L_pk St_L_Hs(n,2)+40+Lidx(L_id)];
                        n;

                        prop_data_L =[(0.05*(St_L_Hs(n,2)+40+Lidx(L_id):St_L_To(n,2)))', Gyro_St_LS_r(St_L_Hs(n,2)+40+Lidx(L_id):St_L_To(n,2))];
                        [~,~,L_max_acc,L_id_acc]= createFit(prop_data_L(:,1),prop_data_L(:,2));
                        
                        n;
                        abs(L_max_acc);
                        Gyro_St_LS_r(St_StanceL(n,2)+L_id_acc-1);
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
                        tro_angle_L_y(n,:) = [L_pk St_L_Hs(n,2)-201+L_id];
                end
                clear P1 P2
                clear Gyro_St_RS_x_r
                clear R_pk L_pk Ridx Lidx R_id L_id
                clear n m R_pks L_pks
                clear L_id_acc R_id_acc Lpks Rpks
                clear prop_data_R prop_data_L  Gyro_St_RS_r  Gyro_St_LS_r
                clear L_max_acc R_max_acc
                
                St_R_To
                St_StanceR
                St_L_To
                
                
    % Calculation of Parameters
    Param_St_L_Hs_Pk = St_StanceL - St_L_Hs;
    Param_St_L_Pk_To = St_L_To    - St_StanceL;
    Param_St_L_Hs_To = St_L_To    - St_L_Hs;
    Param_St_R_Hs_Pk = St_StanceR - St_R_Hs;
    Param_St_R_Pk_To = St_R_To    - St_StanceR;
    Param_St_R_Hs_To = St_R_To    - St_R_Hs;
    clear St_StanceL St_StanceR

    
    
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
    clear St_R_To St_L_To Param_St_L_Hs_R_Hs  Param_St_R_Hs_L_Hs
    clear St_L_Hs St_R_Hs 
    Param_St_L_Hs_Pk;
    Param_St_L_Pk_To
    Param_St_L_Hs_To
    
    Result_St_L = [Param_St_L_Hs_Pk, Param_St_L_Pk_To, Param_St_L_Hs_To, xParam_ST_L_Pk_Pk, Param_x_L_Pk_Tr, Param_y_L_Pk_Tr];
    Result_St_R = [Param_St_R_Hs_Pk Param_St_R_Pk_To Param_St_R_Hs_To xParam_ST_R_Pk_Pk Param_x_R_Pk_Tr Param_y_R_Pk_Tr];
    
    clear Param_x_L_Pk_Tr Param_x_R_Pk_Tr Param_y_L_Pk_Tr Param_y_R_Pk_Tr
    clear Param_St_L_Hs_Pk Param_St_L_Hs_To Param_St_L_Pk_To Param_St_R_Hs_Pk Param_St_R_Hs_To Param_St_R_Pk_To
    clear xParam_ST_L_Pk_Pk xParam_ST_R_Pk_Pk
    
    slope_R = -Result_St_R(:,3)./Result_St_R(:,4)*1000     % Average slope of angular velocity for midstance peak to toe-off
    slope_L = -Result_St_L(:,3)./Result_St_L(:,4)*1000
    prop_time_ratio_R = Result_St_R(:,4)./Result_St_R(:,6); % Propultion to stance time ratio
    prop_time_ratio_L = Result_St_L(:,4)./Result_St_L(:,6);
    
    slope_R_angle_x = Result_St_R(:,9)./Result_St_R(:,10)*1000;     % Average slope of angle x-axis
    slope_L_angle_x = Result_St_L(:,9)./Result_St_L(:,10)*1000;
    
    slope_R_angle_y = Result_St_R(:,11)./Result_St_R(:,12)*1000;     % Average slope of angle y-axis
    slope_L_angle_y = Result_St_L(:,11)./Result_St_L(:,12)*1000;
    
    % Ouptut in structure format
    DataOutput(zz).ID = PatientNameData;
    DataOutput(zz).PropulsionDuration_R =  mean(Result_St_R(:,6))/1000;
    DataOutput(zz).PropulsionDuration_L =  mean(Result_St_L(:,6)/1000);
    DataOutput(zz).PropulsionAcceleration_R = mean(slope_R);
    DataOutput(zz).PropulsionAcceleration_L = mean(slope_L);
    DataOutput(zz).MidstanceSpeed_R = mean(-Result_St_R(:,3));
    DataOutput(zz).MidstanceSpeed_L = mean(-Result_St_L(:,3));
    DataOutput(zz).Speednorm_R = mean(peak_abs_xy_R(:,1));
    DataOutput(zz).Speednorm_L = mean(peak_abs_xy_L(:,1));
    DataOutput(zz).Toeoffspeed_R = mean(peak_StanceR_x(:,1));
    DataOutput(zz).Toeoffspeed_L = mean(peak_StanceL_x(:,1));
    DataOutput(zz).PeakSwingSpeed_R = mean(St_R_PSS);
    DataOutput(zz).PeakSwingSpeed_L = mean(St_L_PSS);
    


    x = Result_St_R(:,6);
%    clear Result_St_L Result_St_R slope_R slope_L
    
    
   %{ 
    output_results(zz,:)=[PatientNameDataDble, ...
        mean(St_R_PSS),std(St_R_PSS),std(St_R_PSS)/mean(St_R_PSS),...                                               % Peak Swing Speed outputs (PSS)
        mean(St_L_PSS),std(St_L_PSS),std(St_L_PSS)/mean(St_L_PSS),...
        mean(-Result_St_R(:,3)),std(-Result_St_R(:,3)),std(-Result_St_R(:,3))/mean(-Result_St_R(:,3)),...           % Midstance peak velocity outputs (delta V) (PMS)
        mean(-Result_St_L(:,3)),std(-Result_St_L(:,3)),std(-Result_St_L(:,3))/mean(-Result_St_L(:,3)),...
        mean(slope_R),std(slope_R),std(smmlope_R)/mean(slope_R),...                                                   % Average slope of angular velocity for midstance peak to toe-off outputs
        mean(slope_L),std(slope_L),std(slope_L)/mean(slope_L),...
        mean(prop_time_ratio_R),std(prop_time_ratio_R),std(prop_time_ratio_R)/mean(prop_time_ratio_R),...           % Propultion/stance time ratio
        mean(prop_time_ratio_L),std(prop_time_ratio_L),std(prop_time_ratio_L)/mean(prop_time_ratio_L),...
        mean(Result_St_R(:,6)),std(Result_St_R(:,6)),std(Result_St_R(:,6))/mean(Result_St_R(:,6)), ...              % Absolute Propulsion time (delta T)
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
    %}
end
clear s zz ans
return;
% Save data
    cd(ResultsDir);
    str = date;
    FileName = ['FrailtyandPropulsion_',str]
    writetable(struct2table(DataOutput), [ResultsDir,FileName,'.xlsx'])
    cd(CurDir)  
    clear str
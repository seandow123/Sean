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

    
%% Directory Location
    % Code Directory
    CodeDir = fullfile(pwd);
    % Raw data directory
    RawDataDir  = 'Z:\Projects BCM\H-38994 BLANKET BCM\Studies\Propulsion\Data\Raw Sensor data\Gait\';
    
    % Directory to same matlab and excel resutls
    ResultsDir  = 'Z:\Projects BCM\H-38994 BLANKET BCM\Studies\Propulsion\Results\';
    
    %%  Load Data for Heel Area
    
    Area_Excel_Dir = 'Z:\Projects BCM\H-38994 BLANKET BCM\Studies\Propulsion\Results\converted Data';
%     pressure_Dir = 'Z:\Personnel\Mohsen\OHI Fscan\Excel Data';
    pressure_Dir = 'Z:\Projects BCM\H-38994 BLANKET BCM\Studies\Propulsion\Data\Analyzed\Fscan Healthy';

%     pressure_Dir = 'Z:\Projects BCM\PropulsionandFrailty\Data\Analyzed\FScanHealthy';
    cd(pressure_Dir);
    [file, path,index] = uigetfile({'*.csv'},...
                          'File Selector');
    cd(path)
    [pressure_Whole,~,~] =xlsread(file);
    [file2, path,index] = uigetfile({'*.csv'},...
                          'File Selector');
    [pressure_WholeTest,~,~] =xlsread(file2);
    
%     [file3, path,index] = uigetfile({'*.csv'},...
%                           'File Selector');
%     [pressure_Heel,~,~] =xlsread(file3);
%     [file4, path,index] = uigetfile({'*.csv'},...
%                           'File Selector');
%     [pressure_HeelTest,~,~] =xlsread(file4);

    %     force = Area_Heel(:,5);
%     figure('Position',[10,10,1800,1800]);
%     plot(force);
%     [points ~] = getpts;
%     points = points *2;
%     close;
    
    pressure = pressure_Whole(:,5);
    figure('Position',[10,10,1800,1800]);
    plot(pressure);
    
%     return;
    
    [points ~] = getpts;
    points = points *2;
    close;
    
%     
    pressureTest = pressure_WholeTest(:,5);
    figure('Position',[10,10,1800,1800]);
    plot(pressureTest);
    
%     return;
    
    [pointsTest ~] = getpts;
    pointsTest = pointsTest *2;
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
        
        
        
Points6 = [2839 2933;3078 3147;3302 3367;3719 3804;3952 4011;4180 4248];


Points26 = [4238 4310;4455 4526;4668 4741;4900 4960;5323 5409;5555 5643];
         
%% Loop for Analyzing Each Data Folder
%   The code is modified to account for follows up inside each other 
%   Optiontime (Baseline - etc.)
for zz = 1:1        % length of selected folder for analysis
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
            
            
            
%             Gyro_test_RS = filtfilt([1,-1],[1,-0.995],-SensorData(StartPt:EndPt,4));
%             Gyro_test_LS = filtfilt([1,-1],[1,-0.995],-SensorData(StartPt:EndPt,24));
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
%                 clear Height 
            
            % Test plot
%                     fig = figure(1);
%                     clf(fig)
%                     pressure_resample = resample(pressure,2,1);
%                     h(1)=subplot(2,1,1); hold on; plot(Gyro_St_RS_r);h(1)=subplot(2,1,1); hold on; plot(pressure_resample);
% %                     xlim([1 length(pressure_resample)])
%                     grid on;
%                     
%                     pressureTest_resample = resample(pressureTest,2,1);
%                     h(2)=subplot(2,1,2);hold on; plot(Gyro_St_LS_r);;h(2)=subplot(2,1,2); hold on; plot(pressureTest_resample);
% %                     xlim([1 length(pressureTest_resample)])
%                     grid on;
%                     
%                     linkaxes(h,'xy')
%                     clear fig h  
            
            
            
            
            
            
%             if(SensorSetup == 2 || SensorSetup ==5)
%                 % Sensor 1 or "Right shin"
%                 Gyro_St_RS(:,1) =  SensorData(StartPt:EndPt,2);	% x for s1
%                 Gyro_St_RS(:,2) =  SensorData(StartPt:EndPt,3);	% y for s1
%                 Gyro_St_RS(:,3) = -SensorData(StartPt:EndPt,4);	% z for s1
%                 
%                 % Sensor 2 or "Left Shin"
%                 Gyro_St_LS(:,1) =  SensorData(StartPt:EndPt,22);	% x for s3
%                 Gyro_St_LS(:,2) =  SensorData(StartPt:EndPt,23);	% y for s3
%                 Gyro_St_LS(:,3) = -SensorData(StartPt:EndPt,24);	% z for s3
%             else
%                 error('Wrong sensor setup')  
%             end

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

        % Test plot
%             x1 = 1;
%             if (x1 == 1)
%                 x2 = max(length(Gyro_St_RS(:,3)));
%             else
%                 x2 = 1000 + x1;
%             end
% 
%             fig = figure(1);
%             clf(fig);
%             h(1)=subplot(2,1,1); plot(Gyro_St_RS(:,3));
%             title([PatientNameData,' Right Shin Gyro about z'])
%             grid on;
%             xlim([x1 x2])
% 
%             h(2)=subplot(2,1,2); plot(Gyro_St_LS(:,3));
%             title([PatientNameData,' Left Shin Gyro about z'])
%             xlim([x1 x2])
%             grid on
%             linkaxes(h,'x')
%             clear Gyro_St_RS Gyro_St_LS x1 x2 fig h  
            
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
                
                % Test plot
%                     fig = figure(1);
%                     h(1)=subplot(2,1,1); hold on; plot(Gyro_St_RS);
%                     ylabel('Angular Velocity')
% 
%                     h(2)=subplot(2,1,2); hold on; plot(Gyro_St_LS);
%                     grid on
%                     linkaxes(h,'x')
%                     ylabel('Angular Velocity')
%                     clear fig h  
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
            
                % Test plot
%                     fig = figure(1);
%                     clf(fig)
%                     h(1)=subplot(2,1,1); plot(Gyro_St_RS_r);
%                     h(1)=subplot(2,1,1); plot(Acc_St_RS_x_r);
%                     h(1)=subplot(2,1,1); plot(Acc_St_RS_y_r);
%                     xlim([1 length(Gyro_St_RS_r)])
%                     grid on;
% 
%                     h(2)=subplot(2,1,2);plot(Gyro_St_LS_r);
%                     xlim([1 length(Gyro_St_RS_r)])
%                     
%                     grid on
%                     linkaxes(h,'xy')
%                     clear fig h  
        
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
                
                Propulsion_point = FindPropulsion(accelometer,[res_St.HsR],[res_St.ToR]);
                Propulsion_point_L = FindPropulsion(accelometer_L,[res_St(1:end-1).HsL],[res_St(2:end).ToL]);
                    % Test plot
                    fig = figure(1);
                    clf(fig)
%                     h(1)=subplot(2,1,1); 
                    hold on; plot(Gyro_St_RS_r);
%                     h(1)=subplot(2,1,1); hold on;
                    plot(accelometer * 400);
                    plot(DG_x * 400,'c');
%                     h(1)=subplot(2,1,1); hold on; plot(DG_x * 400);
%                      h(1)=subplot(2,1,1); hold on; 
                     plot(Propulsion_point,accelometer(Propulsion_point)* 400,'*r','markersize',12);
%                     h(1)=subplot(2,1,1); hold on; plot(-DG_x * 400,'c');
%                      h(1)=subplot(2,1,1); hold on; plot(-DG_y * 400,'k');
                    xlim([1 length(Gyro_St_RS_r)])
                    grid on;

%                     h(2)=subplot(2,1,2);hold on; plot(Gyro_St_LS_r);
%                     h(2)=subplot(2,1,2); hold on; plot(accelometer_L * 400);
% %                     h(1)=subplot(2,1,1); hold on; plot(DG_x * 400);
%                      h(2)=subplot(2,1,2); hold on; plot(Propulsion_point_L,accelometer_L(Propulsion_point_L)* 400,'*r','markersize',12);
%                     xlim([1 length(Gyro_St_RS_r)])
%                     grid on;
                    
%                     linkaxes(h,'xy')
                    clear fig                  
    % Analysis
        Gyro_St_ss = [Gyro_St_RS_r, Gyro_St_LS_r];
        res_St     = GaitAnalyze(Gyro_St_ss, Height);
        clear Height 
        
        
        
        
           % Test plot
           [~, c1] = size(res_St);
           for i = 1:c1-1
               if (file(8) == 'R' || file(7) == 'R')
                   fig = figure(1);
%                    h(1)=subplot(2,1,1); hold on; 
                   plot(res_St(i).HsR,Gyro_St_RS_r(res_St(i).HsR),'.k','markersize',12);
%                    h(1)=subplot(2,1,1); hold on; 
                   plot(res_St(i).ToR,Gyro_St_RS_r(res_St(i).ToR),'*r','markersize',12);
%                    
%                    h(2)=subplot(2,1,2); hold on; plot(res_St(i).HsL,Gyro_St_LS_r(res_St(i).HsL),'.k','markersize',12);
%                    h(2)=subplot(2,1,2); hold on; plot(res_St(i).ToL,Gyro_St_LS_r(res_St(i).ToL),'*r','markersize',12);
                   
                   
%                    pressure_resample_whole = resample(pressure_Whole(:,5),2,1);
%                    pressure_resample_wholeTest = resample(pressure_WholeTest(:,5),2,1);
%                    xlim([1 length(pressurae_resample_whole)])
                   ylim([-300 600]);
%                    linkaxes(h,'xy')
                   
                   
               end
           end
           
           
           pressure_resample_whole = resample(pressure_Whole(:,5),2,1);
                   pressure_resample_wholeTest = resample(pressure_WholeTest(:,5),2,1);

                   
                   % Filter parameters
                Fs = 200;
                filterorder = 11;
                filtercutoff = 15/(Fs);
                filtertype = 'low';
                [b,a] = butter(filterorder,filtercutoff,filtertype);
                clear Fs filterorder filtercutoff filtertype
            % Right shin
            
                pressure_resample_whole = filtfilt(b,a,pressure_resample_whole);
                pressure_resample_wholeTest = filtfilt(b,a,pressure_resample_wholeTest);
           
                pressure_resample_whole = detrend(pressure_resample_whole);
                pressure_resample_wholeTest = detrend(pressure_resample_wholeTest);
           
           pressure_resample_whole = sgolayfilt(pressure_resample_whole,3,21);
            pressure_resample_wholeTest = sgolayfilt(pressure_resample_wholeTest,3,21);
%            if (str2double(file(6:7)) == PatientNameDataDble && str2double(file2(6:7)) == PatientNameDataDble)
%                
%                
%                if (file(8) == 'R')
%                    
%                    if (points > res_St(2).HsR)
%                        h(1)=subplot(2,1,1); hold on; plot(pressure_resample_whole(round(points)-res_St(2).HsR:end),'g');
%                    else
%                        h(1)=subplot(2,1,1); hold on; plot(res_St(2).HsR - round(points):length(pressure_resample_whole)+res_St(2).HsR - round(points)-1,pressure_resample_whole,'g');
%                    end
% %                    
% %                    x = input('number of the peaks: ');
% %                    
% %                    for k = 1:x
% %                        [pointpropul y] = getpts;
% %                        pointofInterest(k,1) = round(pointpropul);
% %                        pointofInterest(k,2) = Gyro_St_RS_r(pointofInterest(k,1));
% %                        h(1)=subplot(2,1,1); hold on; plot(pointofInterest(k,1),pointofInterest(k,2),'Ob','markersize',4);
% %                    end
%                    %                    if (file2(8) == 'L')
%                    if (pointsTest > res_St(1).HsL)
%                        h(2)=subplot(2,1,2); hold on; plot(pressure_resample_wholeTest(round(pointsTest)-res_St(1).HsL:end),'c');
%                    else
%                        h(2)=subplot(2,1,2); hold on; plot(res_St(1).HsL - round(pointsTest):length(pressure_resample_wholeTest)+res_St(1).HsL - round(pointsTest)-1,pressure_resample_wholeTest,'c');
%                    end
%                    
% %                    x = input('number of the peaks: ');
% %                    
% %                    for k = 1:x
% %                        [pointpropul y] = getpts;
% %                        pointofInterest2(k,1) = round(pointpropul);
% %                        pointofInterest2(k,2) = Gyro_St_LS_r(pointofInterest2(k,1));
% %                        h(2)=subplot(2,1,2); hold on; plot(pointofInterest2(k,1),pointofInterest2(k,2),'Ob','markersize',4);
% %                    end
%                    
%                    
%                    %                    else
%                    %                        error('Error. \n Data did not chose correctly.')
%                end
%               
%            end
%                    
%                
%                
           time = 0:.005:(length(Gyro_St_LS_r)-1)/200;
%            timepressure = 0:.005:(length(pressure_resample_Heel)-1)/200;
           time = time';
%            timepressure = timepressure';
           IntegGyro_St_RS_r = cumtrapz(Gyro_St_RS_r);
           IntegGyro_St_LS_r = cumtrapz(Gyro_St_LS_r);
           dGyro_St_RS_rdx = diff(Gyro_St_RS_r(:))./diff(time(:));
           dGyro_St_LS_rdx = diff(Gyro_St_LS_r(:))./diff(time(:));
           d2Gyro_St_RS_rdx = diff(dGyro_St_RS_rdx(:))./diff(time(1:end-1));
           d2Gyro_St_LS_rdx = diff(dGyro_St_LS_rdx(:))./diff(time(1:end-1));
%            dpressureHeeldx = diff(pressure_resample_Heel(:))./diff(timepressure(:));
%            dpressureHeel2dx = diff(pressure_resample_HeelTest(:))./diff(timepressure(:));
%            dGyro_St_RS_rdx = diff([eps; Gyro_St_RS_r(:)])./diff([eps; time(:)]);
%            dGyro_St_LS_rdx = diff([eps; Gyro_St_LS_r(:)])./diff([eps; time(:)]);
%            if (str2double(file(6:7)) == PatientNameDataDble && str2double(file2(6:7)) == PatientNameDataDble)
               
               
                    if (file(8) == 'R' || file(7) == 'R')
                        %%

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
               %%
               
%                    if (pointsTest > res_St(1).ToR)
%                        h(1)=subplot(2,1,1); hold on; plot(pressure_resample_wholeTest(round(pointsTest)-res_St(1).ToR:end),'c');
%                    else
%                        h(1)=subplot(2,1,1); hold on; plot(res_St(1).ToR - round(pointsTest):length(pressure_resample_wholeTest)+res_St(1).ToR - round(pointsTest)-1,pressure_resample_wholeTest,'c');
%                    end
                    
%                     if(PatientNameDataDble ==6)
%                     for l = 1: 6
%                            [VT1best,VT1m1,VT1m2, VT1b1, VT1b2, kbest,kbestMSE,RSSmin, MSEmin, R2best,MSE] =  VT_Vslope_backward_gui(0:(Points6(l,2)-5-Points6(l,1)),Gyro_St_RS_r(Points6(l,1):Points6(l,2)-5),2);
%                            PointforPlot6(l,1) = Points6(l,1) + kbest;
%                            PointforPlot6(l,2) = Gyro_St_RS_r(Points6(l,1) + kbest);
%                     end 
%                     h(1)=subplot(2,1,1); hold on; plot(PointforPlot6(:,1),PointforPlot6(:,2),'*b','markersize',4);
%                    end
% 
%                     if(PatientNameDataDble ==26)
%                     for l = 1: 6
%                            [VT1best,VT1m1,VT1m2, VT1b1, VT1b2, kbest,kbestMSE,RSSmin, MSEmin, R2best,MSE] =  VT_Vslope_backward_gui(0:(Points26(l,2)-5-Points26(l,1)),Gyro_St_RS_r(Points26(l,1):Points26(l,2)-5),2);
%                            PointforPlot26(l,1) = Points26(l,1) + kbest;
%                            PointforPlot26(l,2) = Gyro_St_RS_r(Points26(l,1) + kbest);
%                     end 
%                     h(1)=subplot(2,1,1); hold on; plot(PointforPlot26(:,1),PointforPlot26(:,2),'*b','markersize',4);
%                     end


%%
% select the propulsion based on the first derivitive of the pressure
%                     for u = 2:length(res_St)-1
%                         if (points > res_St(3).HsR)
%                             
%                             
%                         end
%                         
%                         
%                         
%                     end

                   %%
%                    h(1)=subplot(2,1,1); hold on;plot(1:length(Gyro_St_RS_r)-1,diff(Gyro_St_RS_r) * 10,'b');
                    for l = 2:length(res_St)-1
% peak based on the peak and toe 
%                        peakRightStride = max( Gyro_St_RS_r(res_St(l).HsR + 40:res_St(l).ToR));
%                        pointPeakRightStride = find(Gyro_St_RS_r(res_St(l).HsR + 40:res_St(l).ToR) == peakRightStride) + res_St(l).HsR + 30 ;
%                        
%                        h(1)=subplot(2,1,1); hold on;plot(pointPeakRightStride:(res_St(l).ToR-6),diff(Gyro_St_RS_r(pointPeakRightStride:res_St(l).ToR-5)) * 20,'k');
%                        
% %                        peak = findpeaks(diff(Gyro_St_RS_r(pointPeakRightStride:res_St(l).ToR-5)));
%                        [peak,locs,widths,proms] = findpeaks(diff(Gyro_St_RS_r(pointPeakRightStride:res_St(l).ToR-5)));
%                        widths
%                        biggestPeak = find(widths == max(widths))
%                        PointofPeakGyr = find(diff(Gyro_St_RS_r(pointPeakRightStride:res_St(l).ToR-5)) == peak(biggestPeak));
%                        h(1)=subplot(2,1,1); hold on; plot(PointofPeakGyr + pointPeakRightStride,Gyro_St_RS_r(PointofPeakGyr + pointPeakRightStride),'*r','markersize',12);


% peak based on the swing movement in the other side and Heel strike 
                        PointforSwingadjR(l-1,1) = find(Gyro_St_LS_r == res_St(l).PeakSwingSpeedL);
                        PointForHeelStrikeCounterlaterla = res_St(l).HsL;
%                         PointforSwingadjR(l-1,2) = find(Gyro_St_RS_r(res_St(l).HsR:res_St(l).ToR) == PointforSwingadjR(l-1,1)) + res_St(l).HsR ;
%                         h(1)=subplot(2,1,1); hold on;plot(PointforSwingadjR(l-1,1):PointForHeelStrikeCounterlaterla-1,diff(Gyro_St_RS_r(PointforSwingadjR(l-1,1):PointForHeelStrikeCounterlaterla)) * 20,'k');
                        plot(PointforSwingadjR(l-1,1),50,'*k','markersize',12);


                        [peak,locs,widths,proms] = findpeaks(diff(Gyro_St_RS_r(PointforSwingadjR(l-1,1):PointForHeelStrikeCounterlaterla)));
                            biggestPeak = find(widths == max(widths))
                            if(size(peak,1) ~= 0 )
                                
                                PointofPeakGyr = find(diff(Gyro_St_RS_r(PointforSwingadjR(l-1,1):PointForHeelStrikeCounterlaterla)) == peak(biggestPeak));
                                
                                
%                                 h(1)=subplot(2,1,1); hold on; plot(PointofPeakGyr + PointforSwingadjR(l-1,1),Gyro_St_RS_r(PointofPeakGyr + PointforSwingadjR(l-1,1)),'*r','markersize',12);
                                %
                            end
                       
%                        [VT1best,VT1m1,VT1m2, VT1b1, VT1b2, kbest,kbestMSE,RSSmin, MSEmin, R2best,MSE] =  VT_Vslope_backward_gui(0:res_St(l).ToR -pointPeakRightStride -5 ,Gyro_St_RS_r(pointPeakRightStride:res_St(l).ToR-5),2);
%                        PointforPlotR(l-1,1) = pointPeakRightStride + kbest;
%                        PointforPlotR(l-1,2) = Gyro_St_RS_r(pointPeakRightStride + kbest);
%                        h(1)=subplot(2,1,1); hold on; plot(PointforPlotR(l-1,1),PointforPlotR(l-1,2),'*b','markersize',10);
                    end

                    %%
%                     for l = 2:length(res_St)-1
%                        PointforSwingadjR(l-1,1) = find(Gyro_St_LS_r == res_St(l).PeakSwingSpeedL);
%                        PointofToeAdjR(l-1,1) = res_St(l).HsL;
%                        
% %                        pointsStride = find(
%                        
% %                        [pointpropul y] = getpts;
% %                        pointofInterestEnd(k,1) = round(pointpropul);
% %                        pointofInterestEnd(k,2) = Gyro_St_RS_r(pointofInterestEnd(k,1));
%                        
% %                        PointforSwingadjR(l-1,2) = Gyro_St_RS_r(PointforSwingadjR(l-1,1));
% %                        PointofToeAdjR(l-1,2) = Gyro_St_RS_r(PointofToeAdjR(l-1,1));
%                        
% %                        h(1)=subplot(2,1,1); hold on; plot(pointofInterestEnd(l-1,1),pointofInterestEnd(l-1,2),'*k','markersize',10);
%                        
% %                          [~, ~, max_acc,id_acc] = createFit(0:res_St(l).ToR -PointforSwingadjR(l-1,1) -15 ,Gyro_St_RS_r(PointforSwingadjR(l-1,1):res_St(l).ToR -15));
% %                        [VT1best,VT1m1,VT1m2, VT1b1, VT1b2, kbest,kbestMSE,RSSmin, MSEmin, R2best,MSE] =  VT_Vslope_backward_gui(0:res_St(l).ToR -PointforSwingadjR(l-1,1) -15 ,Gyro_St_RS_r(PointforSwingadjR(l-1,1):res_St(l).ToR -15),2);
% %                        [VT1best,VT1m1,VT1m2, VT1b1, VT1b2, kbest,kbestMSE,RSSmin, MSEmin, R2best,MSE] =  VT_Vslope_backward_gui(0:PointofToeAdjR(l-1,1) -PointforSwingadjR(l-1,1) ,Gyro_St_RS_r(PointforSwingadjR(l-1,1):PointofToeAdjR(l-1,1)),2);
%                         
% %                        MSEmin;
% 
%                          h(1)=subplot(2,1,1); hold on;plot(PointforSwingadjR(l-1,1):(res_St(l).ToR-6),diff(Gyro_St_RS_r(PointforSwingadjR(l-1,1):res_St(l).ToR-5)) * 20,'k');
% 
% 
% %                        PointforPlotR(l-1,1) = PointforSwingadjR(l-1,1) + id_acc;
% %                        PointforPlotR(l-1,2) = Gyro_St_RS_r(PointforSwingadjR(l-1,1) + id_acc);
%                      
%                        
% %                        PointforPlotR(l-1,1) = PointforSwingadjR(l-1,1) + kbestMSE;
% %                        PointforPlotR(l-1,2) = Gyro_St_RS_r(PointforSwingadjR(l-1,1) + kbestMSE);
% %                        h(1)=subplot(2,1,1); hold on; plot(PointforPlotR(l-1,1),PointforPlotR(l-1,2),'*b','markersize',10);
%                        
%                        
%                        
%                     end

%                    
%                    x = input('number of the peaks: ');
                    x = length(res_St) - 1;
                   
                   for k = 1:x
                       [pointpropul y] = getpts;
                       pointofInterest(k,1) = round(pointpropul);
%                        pointofInterest(k,1) = round(pointpropul);
%                        pointofInterest(k,2) = Gyro_St_RS_r(pointofInterest(k,1));
%                        [pointpropul y] = getpts;
%                        pointofInterestGyro(k,1) = round(pointpropul);
%                        pointofInterestGyro(k,2) = Gyro_St_RS_r(pointofInterestGyro(k,1));
%                        h(1)=subplot(2,1,1); hold on; plot(pointofInterest(k,1),pointofInterest(k,2),'Ob','markersize',10);
%                        h(1)=subplot(2,1,1); hold on; plot(pointofInterestGyro(k,1),pointofInterestGyro(k,2),'.k','markersize',10);
                   end
                   for u = 1:x
                       pointToR(u,1) = res_St(u).ToR;
                       pointHsR(u,1) = res_St(u).HsR;
                   end
                   PointRight = [pointHsR Propulsion_point(1,1:end-1)' pointofInterest pointToR]
                   
                   %                    if (file2(8) == 'L')
                   %% the plot for the left gait 
%                    if (pointsTest > res_St(2).HsL)
%                        h(2)=subplot(2,1,2); hold on; plot(pressure_resample_wholeTest(round(pointsTest)-res_St(2).HsL:end) * 4,'c');
% %                        h(2)=subplot(2,1,2); hold on; plot(pressure_resample_HeelTest(round(pointsTest)-res_St(2).HsL:end),'c');
% %                        h(2)=subplot(2,1,2); hold on; plot(IntegGyro_St_LS_r,'k');
% %                        h(2)=subplot(2,1,2); hold on; plot(Gyro_St_LS_x_r,'k');
% %                        h(2)=subplot(2,1,2); hold on; plot(dpressureHeel2dx,'c');
%                    else
%                        h(2)=subplot(2,1,2); hold on; plot(res_St(2).HsL - round(pointsTest):length(pressure_resample_wholeTest)+res_St(2).HsL - round(pointsTest)-1,pressure_resample_wholeTest,'c');
% %                        h(2)=subplot(2,1,2); hold on; plot(res_St(2).HsL - round(pointsTest):length(pressure_resample_HeelTest)+res_St(2).HsL - round(pointsTest)-1,pressure_resample_HeelTest,'c');
% %                        h(2)=subplot(2,1,2); hold on; plot(IntegGyro_St_LS_r,'k');
% %                        h(2)=subplot(2,1,2); hold on; plot(Gyro_St_LS_x_r,'k');
% %                        h(2)=subplot(2,1,2); hold on; plot(dpressureHeel2dx,'c');
%                    end
                   
                   %% Using the swing point in the conterlateral side (it did not work)
% %                    for l = 2:length(res_St)-1
% %                        peakRightStride = max( Gyro_St_LS_r(res_St(l).HsL:res_St(l+1).ToL));
% %                        pointPeakRightStride = find(Gyro_St_LS_r(res_St(l).HsL:res_St(l+1).ToL) == peakRightStride) + res_St(l).HsL ;
% %                        [VT1best,VT1m1,VT1m2, VT1b1, VT1b2, kbest,kbestMSE,RSSmin, MSEmin, R2best,MSE] =  VT_Vslope_backward_gui(0:res_St(l+1).ToL -pointPeakRightStride -5 ,Gyro_St_LS_r(pointPeakRightStride:res_St(l+1).ToL-5),2);
% %                        PointforPlotL(l-1,1) = pointPeakRightStride + kbest;
% %                        PointforPlotL(l-1,2) = Gyro_St_LS_r(pointPeakRightStride + kbest);
% % %                        h(2)=subplot(2,1,2); hold on; plot(PointforPlotL(l-1,1),PointforPlotL(l-1,2),'*b','markersize',10);
% %                     end
                   
                   header={'Right alg point','Right alg value','Right Manual point','Right Manual value','Right Toe-off','Left alg point','Left alg value','Left Manual point','Left Manual value','Left Toe-off'};
                   
                   header = {'Heel Strike' , 'propulsion point Acc.' , 'propulsion point Pressure' , 'Toe off'};
                   xlswrite(listDir(s).name,PointRight,'Sheet1','B2');
                   xlswrite(listDir(s).name,header,'Sheet1','B1');
                   
        %% finding the point of propulsion phase based on the pressure graph in the left side
% %                     x = length(res_St) - 2;
% %                    
% %                    for k = 1:x
% %                        [pointpropul y] = getpts;
% %                        pointofInterest2(k,1) = round(pointpropul);
% %                        pointofInterest2(k,2) = Gyro_St_LS_r(pointofInterest2(k,1));
% % %                        h(2)=subplot(2,1,2); hold on; plot(pointofInterest2(k,1),pointofInterest2(k,2),'Ob','markersize',4);
% %                    end
% %                    for u = 1:x
% %                        pointToL(u,1) = res_St(u+2).ToL;
% %                    end
% %                    PointLeft = [PointforPlotL pointofInterest2 pointToL]
% %                    
% %                    pointTotal = [PointRight PointLeft];

               end
              
%            end
           
%            xlswrite(PatientNameData,pointofInterest)
               
               
               
               
               
               
%                
%            if (file(8) == 'R')
%            
%            if (points > res_St(2).HsR)
%                    h(1)=subplot(2,1,1); hold on; plot(pressure_resample_whole(round(points)-res_St(2).HsR:end),'g');
%                else
%                    h(1)=subplot(2,1,1); hold on; plot(res_St(2).HsR - round(points):length(pressure_resample_whole)+res_St(2).HsR - round(points)-1,pressure_resample_whole,'g');
%            end
%                
% %            if (pointsTest > res_St(2).HsR)
% %                    h(1)=subplot(2,1,1); hold on; plot(pressure_resample_wholeTest(round(pointsTest)-res_St(2).HsR:end),'c');
% %                else
% %                    h(1)=subplot(2,1,1); hold on; plot(res_St(2).HsR - round(pointsTest):length(pressure_resample_wholeTest)+res_St(2).HsR - round(pointsTest)-1,pressure_resample_wholeTest,'c');
% %                end
%                
%                x = input('number of the peaks');
%                
%                for k = 1:x
%                   [pointpropul y] = getpts;
%                   pointofInterest(k,1) = round(pointpropul);
%                   pointofInterest(k,2) = Gyro_St_RS_r(pointofInterest(k,1));
%                    h(1)=subplot(2,1,1); hold on; plot(pointofInterest(k,1),pointofInterest(k,2),'Ob','markersize',12);
%                end
%            
%            elseif (file(8) == 'L')
%                if (points > res_St(1).HsL)
%                    h(2)=subplot(2,1,2); hold on; plot(pressure_resample_whole(round(points)-res_St(1).HsL:end),'g');
%                else
%                    h(2)=subplot(2,1,2); hold on; plot(res_St(1).HsL - round(points):length(pressure_resample_whole)+res_St(1).HsL - round(points)-1,pressure_resample_whole,'g');
%                end
%                
% %                if (pointsTest > res_St(2).HsL)
% %                    h(2)=subplot(2,1,2); hold on; plot(pressure_resample_wholeTest(round(pointsTest)-res_St(1).HsL:end),'c');
% %                else
% %                    h(2)=subplot(2,1,2); hold on; plot(res_St(1).HsL - round(pointsTest):length(pressure_resample_wholeTest)+res_St(1).HsL - round(pointsTest)-1,pressure_resample_wholeTest,'c');
% %                end
%            end
           
%            
%            
%            for i =1:c1
%                     fig = figure(1);
%                     h(1)=subplot(2,1,1); hold on; plot(res_St(i).HsR,Gyro_St_RS_r(res_St(i).HsR),'.k','markersize',12);
%                     h(1)=subplot(2,1,1); hold on; plot(res_St(i).ToR,Gyro_St_RS_r(res_St(i).ToR),'*r','markersize',12);
% 
% %                     xlim([1 length(Gyro_St_RS_r)])
% %                     grid on;
% 
%                     h(2)=subplot(2,1,2); hold on; plot(res_St(i).HsL,Gyro_St_LS_r(res_St(i).HsL),'.k','markersize',12);
%                     h(2)=subplot(2,1,2); hold on; plot(res_St(i).ToL,Gyro_St_LS_r(res_St(i).ToL),'*r','markersize',12);
% %                     area_resample_heel = resample(Area_Heel(:,5),2,1) * 10;
%                     pressure_resample_heel = resample(pressure_Heel(:,5),2,1);
%                     pressure_resample_whole = resample(pressure_Whole(:,5),2,1);
%                     xlim([1 length(pressure_resample_heel)])
%                     ylim([-300 700])
% %                     grid on;
%                     
%                     linkaxes(h,'xy')
% %                     clear fig h      
%            end
% %            if (points > res_St(2).HsR)
% %                h(1)=subplot(2,1,1); hold on; plot(area_resample_heel(round(points)-res_St(2).HsR:end));  
% %            else
% %                h(1)=subplot(2,1,1); hold on; plot(res_St(2).HsR - round(points):length(area_resample_heel)+res_St(2).HsR - round(points)-1,area_resample_heel);
% %            end
% 
%             if (points > res_St(2).HsR)
%                h(1)=subplot(2,1,1); hold on; plot(pressure_resample_heel(round(points)-res_St(2).HsR:end));  
%                h(1)=subplot(2,1,1); hold on; plot(pressure_resample_whole(round(points)-res_St(2).HsR:end)); 
%            else
%                h(1)=subplot(2,1,1); hold on; plot(res_St(2).HsR - round(points):length(pressure_resample_heel)+res_St(2).HsR - round(points)-1,pressure_resample_heel);
%                h(1)=subplot(2,1,1); hold on; plot(res_St(2).HsR - round(points):length(pressure_resample_whole)+res_St(2).HsR - round(points)-1,pressure_resample_whole);
%             end
           
            
end
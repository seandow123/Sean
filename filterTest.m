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
        %% Directory Location
    % Code Directory
    CodeDir = fullfile(pwd);
    % Raw data directory
    RawDataDir  = 'Z:\Projects BCM\PropulsionandFrailty\Data\Raw Sensor\Gait\';
    
    
    %% Load Data 
    % Load Data File
        % Location of the data file
        cd(RawDataDir);                                              
        directory = [pwd,filesep];
    % Opend dialog box to select the folder to analyze.
    %   Folder could be analyzied individually or by batch (+1)
        listDir = dir(['*']);    % List of Patient Data (starting with PAD*)                                                               
        Height =160;
    % Select Data Folder(s) for Analysis
        % s: selection - vector of indices of the selected strings (lenght 1
        % in single selection mode. Will be empty([]) when OK is 0.  OK is 1 if you
        % push the OK button, and 0 if you push the cancel button in the
        % figure
        [s,ok] = listdlg('PromptString','Select a folder:', 'ListString',{listDir.name});
        clear ok
        
        pointsIntersst = [2219 2422 2623 2836 3043 3242 3440 3643];
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
%         PatientNameIndex = find(PtsRange(:,1)==PatientNameDataDble);
        clear BegPts len LenNameMax LenName 
        
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

                % Sensor 2 or "Left Shin"
                Gyro_St_LS(:,1) =  SensorData(:,22);	% x for s3
                Gyro_St_LS(:,2) =  SensorData(:,23);	% y for s3
                Gyro_St_LS(:,3) = -SensorData(:,24);	% z for s3
            else
                error('Wrong sensor setup')
            end
% 
%         % Test plot
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
            % Left shin
                Gyro_St_LS_x =  filtfilt([1,-1],[1,-0.995],-SensorData(:,22));
                Gyro_St_LS_y =  filtfilt([1,-1],[1,-0.995],-SensorData(:,23));
                Gyro_St_LS   =  filtfilt([1,-1],[1,-0.995],-SensorData(:,24));  % z
                clear StartPt EndPt
                
%                 % Test plot
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

            Gyro_St_LS_r   = resample(Gyro_St_LS,200,samp_rate);
            Gyro_St_LS_x_r = resample(Gyro_St_LS_x,200,samp_rate);
            Gyro_St_LS_y_r = resample(Gyro_St_LS_y,200,samp_rate);
            clear Gyro_St_RS Gyro_St_RS_x Gyro_St_RS_y
            clear Gyro_St_LS Gyro_St_LS_x Gyro_St_LS_y
            
%                 Test plot
%                     fig = figure(1);
%                     clf(fig)
%                     h(1)=subplot(2,1,1); plot(Gyro_St_RS_r);
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
            % Left shin    
                Gyro_St_LS_x_r = filtfilt(b,a,Gyro_St_LS_x_r);
                Gyro_St_RS_y_r = filtfilt(b,a,Gyro_St_RS_y_r);
                Gyro_St_LS_y_r = filtfilt(b,a,Gyro_St_LS_y_r);
                clear a b
                
                    % Test plot
                    fig = figure(1);
                    clf(fig)
                    plot(Gyro_St_RS_r(2000:3670));
                    xlim([1 length(Gyro_St_RS_r(2000:3670))])
                    grid on;

%                     h(2)=subplot(2,1,2);hold on; plot(Gyro_St_LS_r(2000:36700));
%                     xlim([1 length(Gyro_St_RS_r(2414:3636))])
%                     grid on;
%                     
%                     linkaxes(h,'xy')
                    clear fig h                  
    % Analysis
        Gyro_St_ss = [Gyro_St_RS_r, Gyro_St_LS_r];
        res_St     = GaitAnalyze(Gyro_St_ss, Height);
        clear Height 
        
        
        
       
        
        
        end
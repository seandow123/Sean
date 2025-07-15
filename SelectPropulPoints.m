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
    ExcelDir = 'Z:\Projects BCM\PropulsionandFrailty\Results\';
    cd(ExcelDir)
    [PtsRange,~,~] =xlsread('OHi_Points.xlsx');
    
%% Directory Location
    % Code Directory
    CodeDir = fullfile(pwd);
    % Raw data directory
    RawDataDir  = 'Z:\Projects BCM\PropulsionandFrailty\Data\Raw Sensor\Gait\';
    % Directory to same matlab and excel resutls
    ResultsDir  = 'Z:\Projects BCM\PropulsionandFrailty\Results\';
    
    %%  Load Data for Heel Area
    
    Area_Excel_Dir = 'Z:\Projects BCM\PropulsionandFrailty\Results\converted Data';
    pressure_Dir = 'Z:\Personnel\Mohsen\OHI Fscan\Excel Data';
    cd(pressure_Dir);
    [file, path,index] = uigetfile({'*.csv'},...
                          'File Selector');
    [pressure_Whole,~,~] =xlsread(file);
    [file2, path,index] = uigetfile({'*.csv'},...
                          'File Selector');
    [pressure_WholeTest,~,~] =xlsread(file2);

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
            
            
            
            Gyro_test_RS = filtfilt([1,-1],[1,-0.995],-SensorData(StartPt:EndPt,4));
            Gyro_test_LS = filtfilt([1,-1],[1,-0.995],-SensorData(StartPt:EndPt,24));
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
                clear Height 
            
            % Test plot
                    fig = figure(1);
                    clf(fig)
                    pressure_resample = resample(pressure,2,1);
                    h(1)=subplot(2,1,1); hold on; plot(Gyro_St_RS_r);h(1)=subplot(2,1,1); hold on; plot(pressure_resample);
                    xlim([1 length(pressure_resample)])
                    grid on;
                    
                    pressureTest_resample = resample(pressureTest,2,1);
                    h(2)=subplot(2,1,2);hold on; plot(Gyro_St_LS_r);;h(2)=subplot(2,1,2); hold on; plot(pressureTest_resample);
                    xlim([1 length(pressureTest_resample)])
                    grid on;
                    
                    linkaxes(h,'xy')
                    clear fig h  
            
            
%             
%             
%             
%             
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
%             
%         % Filter data
%             % Right shin
%                 Gyro_St_RS_x =  filtfilt([1,-1],[1,-0.995],-SensorData(StartPt:EndPt,2));
%                 Gyro_St_RS_y =  filtfilt([1,-1],[1,-0.995],-SensorData(StartPt:EndPt,3));
%                 Gyro_St_RS   =  filtfilt([1,-1],[1,-0.995],-SensorData(StartPt:EndPt,4));	% z
%             % Left shin
%                 Gyro_St_LS_x =  filtfilt([1,-1],[1,-0.995],-SensorData(StartPt:EndPt,22));
%                 Gyro_St_LS_y =  filtfilt([1,-1],[1,-0.995],-SensorData(StartPt:EndPt,23));
%                 Gyro_St_LS   =  filtfilt([1,-1],[1,-0.995],-SensorData(StartPt:EndPt,24));  % z
%                 clear StartPt EndPt
%                 
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
%         % Resampling
%             Gyro_St_RS_r   = resample(Gyro_St_RS,200,samp_rate);
%             Gyro_St_RS_x_r = resample(Gyro_St_RS_x,200,samp_rate);
%             Gyro_St_RS_y_r = resample(Gyro_St_RS_y,200,samp_rate);
% 
%             Gyro_St_LS_r   = resample(Gyro_St_LS,200,samp_rate);
%             Gyro_St_LS_x_r = resample(Gyro_St_LS_x,200,samp_rate);
%             Gyro_St_LS_y_r = resample(Gyro_St_LS_y,200,samp_rate);
%             clear Gyro_St_RS Gyro_St_RS_x Gyro_St_RS_y
%             clear Gyro_St_LS Gyro_St_LS_x Gyro_St_LS_y
%             
%                 % Test plot
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
%         
%         % Low Pass filtering
%             % Filter parameters
%                 Fs = 200;
%                 filterorder = 7;
%                 filtercutoff = 15/(Fs/2);
%                 filtertype = 'low';
%                 [b,a] = butter(filterorder,filtercutoff,filtertype);
%                 clear Fs filterorder filtercutoff filtertype
%             % Right shin
%                 Gyro_St_RS_r = filtfilt(b,a,Gyro_St_RS_r);
%                 Gyro_St_LS_r = filtfilt(b,a,Gyro_St_LS_r);
%                 Gyro_St_RS_x_r = filtfilt(b,a,Gyro_St_RS_x_r);
%             % Left shin    
%                 Gyro_St_LS_x_r = filtfilt(b,a,Gyro_St_LS_x_r);
%                 Gyro_St_RS_y_r = filtfilt(b,a,Gyro_St_RS_y_r);
%                 Gyro_St_LS_y_r = filtfilt(b,a,Gyro_St_LS_y_r);
%                 clear a b
%                 
%                     % Test plot
%                     fig = figure(1);
%                     clf(fig)
%                     h(1)=subplot(2,1,1); hold on; plot(Gyro_St_RS_r);
%                     xlim([1 length(Gyro_St_RS_r)])
%                     grid on;
% 
%                     h(2)=subplot(2,1,2);hold on; plot(Gyro_St_LS_r);
%                     xlim([1 length(Gyro_St_RS_r)])
%                     grid on;
%                     
%                     linkaxes(h,'xy')
%                     clear fig h                  
%     % Analysis
%         Gyro_St_ss = [Gyro_St_RS_r, Gyro_St_LS_r];
%         res_St     = GaitAnalyze(Gyro_St_ss, Height);
%         clear Height 
%         
%         
%         
%         
%            % Test plot
%            [~, c1] = size(res_St);
%            for i = 2:c1
%                if (file(8) == 'R')
%                    fig = figure(1);
%                    h(1)=subplot(2,1,1); hold on; plot(res_St(i).HsR,Gyro_St_RS_r(res_St(i).HsR),'.k','markersize',12);
%                    h(1)=subplot(2,1,1); hold on; plot(res_St(i).ToR,Gyro_St_RS_r(res_St(i).ToR),'*r','markersize',12);
%                    pressure_resample_whole = resample(pressure_Whole(:,5),2,1);
%                    pressure_resample_wholeTest = resample(pressure_WholeTest(:,5),2,1);
%                    xlim([1 length(pressure_resample_whole)])
%                    ylim([-300 700]);
%                    linkaxes(h,'xy')
%                    
%                    
%                    %            elseif (file(8) == 'L')
%                    %                fig = figure(1);
%                    h(2)=subplot(2,1,2); hold on; plot(res_St(i).HsL,Gyro_St_LS_r(res_St(i).HsL),'.k','markersize',12);
%                    h(2)=subplot(2,1,2); hold on; plot(res_St(i).ToL,Gyro_St_LS_r(res_St(i).ToL),'*r','markersize',12);
%                    pressure_resample_whole = resample(pressure_Whole(:,5),2,1);
%                    pressure_resample_wholeTest = resample(pressure_WholeTest(:,5),2,1);
%                    xlim([1 length(pressure_resample_whole)])
%                    ylim([-300 700]);
%                    linkaxes(h,'xy')
%                    
%                    
%                end
%            end
%            
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
%                    
%                    x = input('number of the peaks: ');
%                    
%                    for k = 1:x
%                        [pointpropul y] = getpts;
%                        pointofInterest(k,1) = round(pointpropul);
%                        pointofInterest(k,2) = Gyro_St_RS_r(pointofInterest(k,1));
%                        h(1)=subplot(2,1,1); hold on; plot(pointofInterest(k,1),pointofInterest(k,2),'Ob','markersize',4);
%                    end
%                    %                    if (file2(8) == 'L')
%                    if (pointsTest > res_St(2).HsL)
%                        h(2)=subplot(2,1,2); hold on; plot(pressure_resample_wholeTest(round(pointsTest)-res_St(2).HsL:end),'c');
%                    else
%                        h(2)=subplot(2,1,2); hold on; plot(res_St(2).HsL - round(pointsTest):length(pressure_resample_wholeTest)+res_St(2).HsL - round(pointsTest)-1,pressure_resample_wholeTest,'c');
%                    end
%                    
%                    x = input('number of the peaks: ');
%                    
%                    for k = 1:x
%                        [pointpropul y] = getpts;
%                        pointofInterest2(k,1) = round(pointpropul);
%                        pointofInterest2(k,2) = Gyro_St_LS_r(pointofInterest2(k,1));
%                        h(2)=subplot(2,1,2); hold on; plot(pointofInterest2(k,1),pointofInterest2(k,2),'Ob','markersize',4);
%                    end
%                    
%                    
%                    %                    else
%                    %                        error('Error. \n Data did not chose correctly.')
%                end
%                
%                
%                
%                
%                
%                
%            end
%                    
%                
%                
               
               
               
               
               
               
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
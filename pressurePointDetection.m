%--------------------------------------------------------------------------
% Project Frailty and Propulsion
%   The purpose of the function is to calculate the propulsion data using
%   the legsys sensor
% Author:  Mohsen Zahiri
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
        %% Pressure data
%             pressure_Dir = 'Z:\Personnel\Mohsen\OHI Fscan\Excel Data';

    pressure_Dir = 'Z:\Projects BCM\H-38994 BLANKET BCM\Studies\Fscan Frailty\Data\Raw Sensor data\Fscan Healthy Excel';
    
    %% 
    
    
    
    cd(pressure_Dir);
    [file, path,index] = uigetfile({'*.csv'},...
        'File Selector');
    [pressure_Whole,~,~] =xlsread(file);
    
    pressure_resample_whole = resample(pressure_Whole(:,5),2,1);
  
    % Filter parameters
    Fs = 200;
    filterorder = 11;
    filtercutoff = 15/(Fs);
    filtertype = 'low';
    [b,a] = butter(filterorder,filtercutoff,filtertype);
    clear Fs filterorder filtercutoff filtertype
    % Right shin
    
    pressure_resample_whole = filtfilt(b,a,pressure_resample_whole);
    
    pressure_resample_whole = detrend(pressure_resample_whole);
    
    pressure_resample_whole = sgolayfilt(pressure_resample_whole,3,21);
    
    per10 = prctile(pressure_resample_whole,70);
    per90 = prctile(pressure_resample_whole,90);
    
    
    h(1)=subplot(2,1,1); hold on; plot(pressure_resample_whole * 4,'g'); 
    plot(diff(diff(pressure_resample_whole)) * 200 , 'b')
    plot(per10 * ones(size(pressure_resample_whole)) * 4 , 'b')
    
    
    
    
    
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
%%
% a = [50
% 50
% 50
% 50
% 50
% 50
% 50
% 50
% 50
% 50
% 51
% 50
% 50
% 50
% 51
% 50
% 50
% 50
% 50
% 51
% 51
% 51
% 51
% 51
% 51
% 51
% 51
% 51
% 51
% 51
% 51
% 52
% 52
% 52
% 53
% 53
% 54
% 54
% 55
% 56
% 57
% 58
% 59
% 60
% 61
% 63
% 64
% 66
% 68
% 69
% 69
% 69
% 69
% 69
% 67
% 66
% 65
% 63
% 61
% 59
% 57
% 55
% 53
% 51
% 49
% 46
% 43
% 40
% 37
% 35
% 33
% 33
% 32
% 31
% 30
% 29
% 29
% 27
% 25
% 23
% 22
% 21
% 20
% 21
% 21
% 20
% 20
% 20
% 20
% 20
% 20
% 19
% 19
% 18
% 18
% 18
% 18
% 17
% 17
% 17
% 17
% 17
% 16
% 17
% 16
% 16
% 17
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 18
% 22
% 30
% 42
% 55
% 70
% 81
% 95
% 101
% 104
% 100
% 94
% 91
% 90
% 88
% 87
% 87
% 87
% 84
% 82
% 79
% 78
% 77
% 75
% 73
% 72
% 71
% 68
% 68
% 66
% 66
% 64
% 63
% 62
% 62
% 61
% 61
% 60
% 60
% 60
% 60
% 59
% 59
% 59
% 59
% 60
% 60
% 61
% 63
% 64
% 66
% 69
% 72
% 74
% 77
% 79
% 82
% 83
% 85
% 87
% 88
% 89
% 91
% 93
% 93
% 94
% 93
% 90
% 87
% 82
% 76
% 66
% 57
% 49
% 44
% 37
% 32
% 28
% 24
% 21
% 17
% 17
% 17
% 17
% 16
% 16
% 16
% 16
% 16
% 16
% 17
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 15
% 16
% 16
% 15
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 15
% 15
% 15
% 15
% 15
% 16
% 16
% 16
% 17
% 20
% 29
% 46
% 61
% 76
% 102
% 113
% 115
% 116
% 108
% 103
% 100
% 98
% 99
% 98
% 97
% 91
% 85
% 78
% 74
% 71
% 69
% 66
% 63
% 60
% 59
% 57
% 56
% 56
% 56
% 56
% 56
% 57
% 57
% 58
% 59
% 59
% 61
% 63
% 66
% 68
% 71
% 74
% 77
% 81
% 86
% 87
% 90
% 95
% 95
% 100
% 102
% 103
% 103
% 105
% 108
% 106
% 106
% 102
% 101
% 94
% 87
% 76
% 66
% 56
% 46
% 38
% 32
% 26
% 22
% 19
% 17
% 16
% 17
% 16
% 17
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 17
% 16
% 16
% 16
% 17
% 17
% 17
% 17
% 17
% 17
% 17
% 17
% 18
% 18
% 17
% 16
% 17
% 17
% 17
% 18
% 18
% 17
% 17
% 17
% 16
% 17
% 17
% 18
% 27
% 47
% 68
% 94
% 114
% 119
% 119
% 110
% 102
% 103
% 101
% 98
% 96
% 90
% 84
% 80
% 77
% 74
% 72
% 70
% 66
% 64
% 62
% 60
% 59
% 58
% 57
% 56
% 56
% 56
% 57
% 57
% 57
% 58
% 59
% 61
% 61
% 64
% 67
% 69
% 72
% 76
% 81
% 83
% 88
% 92
% 94
% 97
% 97
% 98
% 100
% 104
% 107
% 108
% 109
% 104
% 103
% 101
% 93
% 83
% 73
% 64
% 55
% 47
% 38
% 32
% 26
% 21
% 18
% 15
% 16
% 15
% 15
% 16
% 15
% 15
% 16
% 16
% 16
% 16
% 16
% 17
% 17
% 17
% 17
% 17
% 17
% 17
% 16
% 17
% 17
% 17
% 16
% 16
% 17
% 17
% 16
% 16
% 16
% 16
% 17
% 17
% 17
% 17
% 17
% 18
% 16
% 19
% 28
% 44
% 63
% 83
% 103
% 115
% 110
% 104
% 95
% 94
% 91
% 89
% 85
% 80
% 78
% 75
% 74
% 72
% 70
% 69
% 66
% 64
% 62
% 60
% 58
% 57
% 56
% 56
% 55
% 55
% 56
% 56
% 57
% 59
% 59
% 59
% 62
% 64
% 66
% 69
% 72
% 77
% 81
% 84
% 88
% 91
% 93
% 96
% 98
% 98
% 102
% 103
% 106
% 104
% 102
% 104
% 104
% 99
% 91
% 82
% 69
% 60
% 50
% 41
% 34
% 29
% 23
% 20
% 17
% 16
% 15
% 15
% 16
% 16
% 16
% 16
% 17
% 16
% 16
% 16
% 16
% 15
% 16
% 16
% 17
% 16
% 16
% 16
% 16
% 16
% 16
% 18
% 17
% 17
% 17
% 17
% 16
% 17
% 17
% 17
% 17
% 17
% 16
% 17
% 18
% 17
% 17
% 17
% 17
% 21
% 32
% 50
% 63
% 82
% 101
% 115
% 123
% 123
% 122
% 110
% 104
% 104
% 105
% 103
% 99
% 94
% 91
% 84
% 78
% 73
% 70
% 67
% 64
% 62
% 60
% 60
% 60
% 59
% 59
% 60
% 62
% 64
% 66
% 68
% 70
% 72
% 74
% 77
% 79
% 81
% 83
% 86
% 88
% 90
% 91
% 94
% 96
% 96
% 96
% 98
% 98
% 98
% 99
% 100
% 100
% 99
% 94
% 87
% 79
% 69
% 61
% 53
% 44
% 36
% 31
% 24
% 19
% 16
% 15
% 15
% 15
% 15
% 15
% 15
% 15
% 15
% 15
% 15
% 16
% 17
% 16
% 17
% 16
% 17
% 16
% 17
% 17
% 16
% 16
% 16
% 16
% 17
% 16
% 16
% 17
% 17
% 16
% 17
% 17
% 18
% 17
% 16
% 16
% 17
% 17
% 19
% 29
% 47
% 71
% 89
% 108
% 111
% 109
% 107
% 97
% 94
% 95
% 96
% 90
% 85
% 81
% 76
% 73
% 70
% 66
% 64
% 60
% 57
% 56
% 54
% 52
% 51
% 50
% 50
% 51
% 50
% 51
% 51
% 51
% 52
% 53
% 54
% 55
% 57
% 59
% 61
% 64
% 67
% 72
% 75
% 79
% 84
% 87
% 91
% 95
% 98
% 99
% 102
% 103
% 103
% 100
% 97
% 91
% 83
% 73
% 64
% 55
% 46
% 37
% 31
% 25
% 21
% 18
% 16
% 15
% 16
% 15
% 15
% 15
% 16
% 16
% 16
% 15
% 16
% 16
% 16
% 15
% 16
% 16
% 15
% 15
% 15
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 17
% 15
% 16
% 16
% 16
% 17
% 17
% 17
% 17
% 16
% 15
% 17
% 23
% 37
% 54
% 75
% 98
% 110
% 109
% 103
% 95
% 93
% 91
% 89
% 84
% 81
% 78
% 76
% 74
% 73
% 70
% 66
% 62
% 60
% 57
% 55
% 54
% 54
% 53
% 53
% 53
% 54
% 55
% 56
% 57
% 59
% 60
% 61
% 62
% 64
% 66
% 68
% 70
% 73
% 78
% 81
% 84
% 87
% 90
% 94
% 97
% 100
% 103
% 105
% 108
% 111
% 108
% 109
% 103
% 99
% 89
% 77
% 67
% 57
% 47
% 39
% 31
% 26
% 21
% 18
% 16
% 16
% 16
% 15
% 15
% 15
% 15
% 15
% 15
% 15
% 15
% 15
% 15
% 15
% 15
% 16
% 15
% 15
% 16
% 15
% 15
% 15
% 15
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 17
% 18
% 19
% 30
% 49
% 71
% 95
% 110
% 118
% 112
% 98
% 94
% 92
% 87
% 86
% 83
% 81
% 80
% 77
% 76
% 73
% 69
% 66
% 63
% 60
% 58
% 56
% 54
% 52
% 51
% 51
% 51
% 51
% 52
% 52
% 53
% 54
% 55
% 57
% 60
% 62
% 63
% 66
% 69
% 75
% 80
% 85
% 92
% 94
% 98
% 100
% 103
% 105
% 108
% 110
% 110
% 109
% 110
% 103
% 95
% 84
% 71
% 62
% 53
% 43
% 35
% 28
% 23
% 19
% 16
% 16
% 16
% 15
% 15
% 15
% 15
% 16
% 15
% 15
% 15
% 15
% 15
% 16
% 16
% 16
% 16
% 15
% 16
% 16
% 16
% 16
% 16
% 15
% 15
% 16
% 15
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 16
% 15
% 16
% 16
% 19
% 30
% 48
% 70
% 87
% 105
% 114
% 106
% 98
% 94
% 92
% 89
% 87
% 82
% 77
% 74
% 71
% 70
% 67
% 64
% 61
% 59
% 57
% 56
% 56
% 56
% 56
% 59
% 60
% 63
% 65
% 68
% 68
% 70
% 71
% 73
% 74
% 76
% 76
% 78
% 81
% 84
% 87
% 90
% 94
% 98
% 103
% 104
% 111
% 110
% 115
% 115
% 118
% 120
% 123
% 121
% 118
% 117
% 110
% 103
% 93
% 79
% 66
% 56
% 46
% 39
% 32
% 25
% 21
% 18
% 16
% 16
% 16
% 16
% 16
% 15
% 16
% 15
% 15
% 15
% 15
% 15
% 15
% 15
% 16
% 16
% 15
% 16
% 15
% 16
% 17
% 17
% 17
% 17
% 16
% 16
% 16
% 16
% 16
% 17
% 16
% 16
% 15
% 17
% 17
% 16
% 22
% 36
% 54
% 76
% 95
% 113
% 110
% 103
% 96
% 95
% 94
% 92
% 88
% 83
% 80
% 77
% 74
% 71
% 69
% 67
% 64
% 62
% 60
% 59
% 57
% 56
% 55
% 54
% 55
% 54
% 54
% 55
% 57
% 59
% 60
% 63
% 66
% 68
% 72
% 73
% 75
% 77
% 78
% 81
% 82
% 85
% 85
% 88
% 89
% 91
% 92
% 95
% 97
% 98
% 99
% 103
% 102
% 104
% 99
% 94
% 92
% 85
% 78
% 67
% 58
% 48
% 42
% 35
% 31
% 27
% 24
% 21
% 20
% 18
% 17
% 17
% 16
% 16
% 15
% 16
% 16
% 16
% 16
% 16
% 15
% 17
% 16
% 17
% 16
% 16
% 16
% 16
% 16
% 16
% 15
% 15
% 16
% 16
% 15
% 16
% 16
% 18
% 18
% 17
% 19
% 19
% 21
% 23
% 25
% 26
% 28
% 29
% 30
% 30
% 31
% 31
% 31
% 32
% 33
% 33
% 34
% 35
% 35
% 36
% 37
% 37
% 38
% 39
% 39
% 40
% 40
% 40
% 40
% 41
% 41
% 41
% 41
% 41
% 41
% 41
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 43
% 42
% 43
% 42
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 44
% 44
% 44
% 44
% 44
% 44
% 44
% 44
% 44
% 44
% 44
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 46
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 46
% 46
% 45
% 45
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 45
% 46
% 45
% 45
% 45
% 45
% 45
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 46
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 46
% 46
% 46
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 44
% 45
% 44
% 44
% 44
% 44
% 44
% 44
% 44
% 44
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 42
% 42
% 43
% 42
% 42
% 42
% 43
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 41
% 42
% 42
% 42
% 42
% 42
% 41
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 43
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 43
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 42
% 42
% 42
% 42
% 43
% 42
% 43
% 43
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 42
% 43
% 43
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 43
% 42
% 42
% 42
% 42
% 43
% 43
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 42
% 43
% 43
% 42
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 42
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 43
% 44
% 44
% 44
% 44
% 44
% 44
% 44
% 44
% 44
% 44
% 44
% 44
% 44
% 44
% 44
% 44
% 44
% 44
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 45
% 46
% 45
% 46
% 46
% 46
% ]
% 
% plot(a)
plot(Gyro_St_LS_r)
%% Directory Location
% Code Directory
CodeDir = fullfile(pwd);
% Raw data directory
% RawDataDir  = 'Z:\Projects BCM\H-38994 BLANKET BCM\Studies\Fscan Frailty\Data\Converted\Gait';

RawDataDir  = 'Z:\Projects BCM\H-38994 BLANKET BCM\Studies\Propulsion\Data\Raw Sensor data\Gait';

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
Gyro_St_ss = [Gyro_St_RS_r, Gyro_St_LS_r]; %left and right sides
res_St     = GaitAnalyze(Gyro_St_ss, Height); %running function code

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
Test_No_wave  = Acc_St_RS_x_r;
Test_No_wave_y  = Acc_St_RS_y_r;
Acc_St_RS_x_r = filtfilt(b,a,Acc_St_RS_x_r);
Acc_St_RS_y_r = filtfilt(b,a,Acc_St_RS_y_r);
[C,L] = wavedec(Acc_St_RS_x_r,5,'db5');
DG_x = wrcoef('a',C,L,'db5',5);


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
[C,L] = wavedec(Acc_St_LS_x_r,5,'db5');
DG_x_L = wrcoef('a',C,L,'db5',5);
[C,L] = wavedec(Acc_St_LS_y_r,4,'db5');
DG_y_L = wrcoef('a',C,L,'db5',4);
accelometer_L = sqrt(DG_x_L .^ 2 + DG_y_L .^2)/2;
clear a b

%% Define the propulsion phase by looking at the accelometer
Propulsion_point_R = FindPropulsionXAccel(DG_x,[res_St.HsR],[res_St.ToR]);
Propulsion_point_L = FindPropulsionXAccel(DG_x_L,[res_St(1:end-1).HsL],[res_St(2:end).ToL]);

%%%%%%%% Legsys process%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot the walking signal and define the Toe off and heel strike


% Fscan Data for healthy subjects
% pressure_Dir = 'Z:\Projects BCM\H-38994 BLANKET BCM\Studies\Fscan Frailty\Data\Raw Sensor data\Fscan Healthy';

% Fscan data for OHI subjects
pressure_Dir = 'Z:\Projects BCM\H-38994 BLANKET BCM\Studies\Propulsion\Data\Analyzed\Fscan OHI';

cd(pressure_Dir);
[file, path,index] = uigetfile({'*.csv'},...
    'File Selector');
% [file_Heel_pressure, path,index] = uigetfile({'*.csv'},...
%     'File Selector');
cd(path)
[pressure_Whole,~,~] =xlsread(file);
% [pressure_Heel,~,~] =xlsread(file_Heel_pressure);

pressure = pressure_Whole(:,5);
figure('Position',[10,10,1800,1800]);
plot(pressure);

[points ~] = getpts;
points = points *2;
close;
pressure_resample_whole = resample(pressure_Whole(:,5),2,1);
% pressure_resample_Heel = resample(pressure_Heel(:,5),2,1);


fig = figure(1);
clf(fig)
hold on; plot(Gyro_St_RS_r);
% plot(accelometer * 400);
% plot(DG_x * 400, 'c');
% plot(DG_x * 400, 'c');

% plot(Test_No_wave * 200, 'k')
% plot(Test_No_wave_y * 100, 'r')
% plot(DG_y * 50, 'k');
plot(Propulsion_point_R,accelometer(Propulsion_point_R)* 400,'*r','markersize',12);
xlim([1 length(Gyro_St_RS_r)])
grid on;
clear fig

% Test plot
[~, c1] = size(res_St);
for i = 1:c1-1
    if (file(8) == 'R' || file(7) == 'R')
        fig = figure(1);
        %                    h(1)=subplot(2,1,1); hold on;
        plot(res_St(i).HsR,Gyro_St_RS_r(res_St(i).HsR),'.k','markersize',12);
        %                    h(1)=subplot(2,1,1); hold on;
        plot(res_St(i).ToR,Gyro_St_RS_r(res_St(i).ToR),'*r','markersize',12);
%         ylim([-150 400]);
    end
end


if (points > res_St(3).HsR)
    
    %                             h(1)=subplot(2,1,1); hold on;
    plot(pressure_resample_whole(round(points)-res_St(3).HsR:end),'g');
%     plot(pressure_resample_Heel(round(points)-res_St(3).HsR:end),'g');

else
    %                             h(1)=subplot(2,1,1); hold on;
    plot(res_St(3).HsR - round(points):length(pressure_resample_whole)+res_St(3).HsR - round(points)-1,pressure_resample_whole,'g');
%      plot(res_St(3).HsR - round(points):length(pressure_resample_Heel)+res_St(3).HsR - round(points)-1,pressure_resample_Heel,'g');
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


%% defining the position of the heel off in the pressure plot


[pointsPre ~] = getpts

for i = 1 : length(pointsPre)
    if (points > res_St(3).HsR)
        plot(round(pointsPre(i)),pressure_resample_whole(round(pointsPre(i))+round(points-res_St(3).HsR)) * 4,'.k','markersize',12);
        result(i).PresR = pointsPre(i);
    else
        plot(round(pointsPre(i)),pressure_resample_whole(round(pointsPre(i))-(round(pointsPre(i))-round(points))) * 4,'.k','markersize',12);
        result(i).PresR = pointsPre(i);
    end
    
end
% prompt = 'Number of Propulsion should be ignored: ';
% 
% num = input(prompt)
% 
% sum = 0;
% for i = num+1:size(res_St,2)
%     ddiff = abs(result(i).PresR-result(i).Pro_R)
%     if ~isempty(ddiff)
%         sum = sum + ddiff;
%     end
% end
% result(num+1).SumErr = sum;


answer = questdlg('Does the plot and points look ok?', ...
    'check the plot', ...
    'No','Yes','Yes');
% Handle response
switch answer
    case 'No'
        
    case 'Yes'
        SavingDir = 'Z:\Projects BCM\H-38994 BLANKET BCM\Studies\Propulsion\Results\PropulsionPoints\PropulsionPint_WholePlantaePressure\Old data-3-12-2019\New Folder';
%         SavingDir = "Z:\Projects BCM\H-38994 BLANKET BCM\Studies\Propulsion\Results\PropulsionPoints\PropulsionPint_WholePlantaePressure\Old data-3-12-2019";
        cd(SavingDir);
        save(strcat(file(1:end-4),'.mat'),'result');
end
    



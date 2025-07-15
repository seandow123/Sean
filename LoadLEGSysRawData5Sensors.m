function [Data_Sensor, samp_rate, dirname] = LoadLEGSysRawDataV3(dirname)
% Mohsen Zahiri and Hung  Nguyen
% Modified code to skip when there is not data
% cd([dirname,filesep])
% dir
%function [quatdata, sensordata, samp_rate] = LoadLEGSysRawData(dirname)
%
% Extract raw data from LEGSys
% returns quaternions and sensor data recorded by LEGSys
% 
% Input:
%   dirname: name of directory containing LEGSys data to load (string)
%            us [] to select directory manually.
%
% Data format:

%  The Legsys software was programed based on  5 sensors. the quaternion data will be extracted from
%  the output raw data and included 20 columns. each 4 column shows the
%  quaternion data for each sensor. if we use less
%  sensor like 2 in PAD study . the 4 first
%  coulmns are for the first sensor and coulmns from 5:8 are zero (The program assume that the second sensot is skiped) and 9:12
%  shows the data for second sensor in the experiment.
%
% quatdata = [q1 q2 q3 q4 q5]
%    each quaternion is 4-tuple
%        [w x y z]
%
% sensor data = [s1 s2 s3 s4 s5]
%    each sensor value is 10-tuple:
%        [accx accy accz magx magy magz gyrox gyroy gyroz temp]
%   s1 = Right Shank
%   s2 = Right Thigh
%   s3 = Left Shank
%   s4 = Left Thigh
%   s5 = Waist
%
% samp_rate = sample rate in Hz
%modifications
%JGwin 02/06/2013 - added optional input path


NUM_SENSORS = 5;           %
CHANNELS_PER_SENSOR = 10;  %acc, mag, gyro, temp
Data_Sensor = [];
% Select the file
if nargin < 1
    dirname = uigetdir;       %select LEGSys raw data directory
else
    listdir = dir([dirname filesep '*.bin']);
    if isempty(listdir)
        disp(' ')
        disp(['No data were found - check folder again'])
        quatdata    = 'error';
        sensordata  = {}; 
        samp_rate   = 0;
        dirname     = 'error';
        return;
%         disp(['This folder appear to be empty, manually select data folder again'])
%         dirname = uigetdir;
%         % check again - and terminate if no data is found
%         listDirManual = dir([dirname filesep '*.bin']);
%         if(isempty(listDirManual))
%             error('No data were available - Check Data again')
%         end
    end
end

% cd(dirname);
fsensor = [dirname filesep 'sensor.bin'];
foutput = [dirname filesep 'output.bin'];

quatdata=import_b(foutput, 'float', 4*NUM_SENSORS);
% quatdata = 'error';
sensordata=import_b(fsensor, 'float', CHANNELS_PER_SENSOR*NUM_SENSORS);

% the data is saves based on 5 sensors, to pick up the data (Acc, Gyro, and Quat.) from two
% sensors, and manage them in one arraye, we used:

Data_Sensor(:,2:4)  = sensordata(:,7:9);
Data_Sensor(:,5:7)  = sensordata(:,1:3);
Data_Sensor(:,8:11) = quatdata(:,1:4);

Data_Sensor(:,12:14) = sensordata(:,17:19);
Data_Sensor(:,15:17) = sensordata(:,11:13);
Data_Sensor(:,18:21) = quatdata(:,5:8);

Data_Sensor(:,22:24) = sensordata(:,27:29);
Data_Sensor(:,25:27) = sensordata(:,21:23);
Data_Sensor(:,28:31) = quatdata(:,9:12);

Data_Sensor(:,32:34) = sensordata(:,37:39);
Data_Sensor(:,35:37) = sensordata(:,31:33);
Data_Sensor(:,38:41) = quatdata(:,13:16);

Data_Sensor(:,42:44) = sensordata(:,47:49);
Data_Sensor(:,45:47) = sensordata(:,41:43);
Data_Sensor(:,48:51) = quatdata(:,17:20);

size(Data_Sensor);
A = [0:.01:size(sensordata,1)/100];
size(A);
Data_Sensor(:,1) = [0:.01:size(sensordata,1)/100-.01];

function [data] = import_b(filename, type, nchannel)
  data = [];
  
  fid = fopen(filename, 'rb');
  buffer = fread(fid,4,'ushort');
  if buffer(1) == 777 && buffer(2) == 0
      %this means that this is the new .bin format
      samp_rate = buffer(3);
      buffer = fread(fid, inf, type, 0, 'ieee-le');
      fclose(fid);
  else
      %this is the old format
      samp_rate = 100;
      fclose(fid);
      fid = fopen(filename, 'rb');
      buffer = fread(fid, inf, type, 0, 'ieee-le');
      fclose(fid);
  end
  
  buffer = buffer(1:end-mod(size(buffer,1),nchannel),1);
  data = reshape(buffer, nchannel , size(buffer,1)/nchannel)';

end
end








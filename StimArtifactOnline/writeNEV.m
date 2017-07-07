function [] = writeNEV(data, packetWidth, filename, mapfileName, comments )
% this function takes spike data and timing as inputs, as well as meta
% information and writes a .NEV file with the given filepath and filename

timeSpikes = data.ts;
waveformSpikes = data.waveforms/0.254; % 1 bit is 254nV 
chanSpikes = data.elec;

% open the NEV file -- add .nev extension if necessary
[pathstr,fname,fext]=fileparts(filename);
fid = fopen(strcat(fname,'.NEV'),'Wb');

%% Section 1 - header basic info.
numExtendedHeaders = 0;
bytesHeader = 0;
headerInfoWrite = {};
precision = {};
% 1. File type ID, 8 bytes, always set as "NEURALEV" for neural events
headerInfoWrite{1,1} = 'NEURALEV';
bytesHeader = bytesHeader + 8;
precision{1,1} = 'char';
% 2. File spec, 2 bytes, major and minor revision numbers -- 0x0203 for 2.3
headerInfoWrite{end+1,1} = 770;
bytesHeader = bytesHeader + 2;
precision{end+1,1} = 'bit16';
% 3. Additional Flags, 2 Bytes, bit 0: set as 1 if all waveforms are 16
%       bit, else if a mixture is expected. All other bits are set to 0
headerInfoWrite{end+1,1} = 1;
bytesHeader = bytesHeader + 2;
precision{end+1,1} = 'bit16';
% 4. Bytes in headers, 4 bytes -- total number of bytes in both headers.
%       Used as a zero idx reference to the first data packet
headerInfoWrite{4} = bytesHeader; % this will updated last
bytesHeader = bytesHeader + 4;
precision{4,1} = 'bit32';
% 5. Bytes in data packets, 4 bytes, length of fixed width data packets in
%       the data section of the file. Must be between 12 and 256, multiple 
%       of 4 (104 will be used for spikes only)
headerInfoWrite{end+1,1} = packetWidth;
bytesHeader = bytesHeader + 4;
precision{end+1,1} = 'bit32';
% 6. Time resolution of time stamps, 4 bytes, frequency of global clock.
%       Will be set to 30,000
headerInfoWrite{end+1,1} = 30000;
bytesHeader = bytesHeader + 4;
precision{end+1,1} = 'bit32';
% 7. Time resolution of samples, 4 bytes, sampling frequency used to
%       digitize neural waveforms. Will be set to 30,000
headerInfoWrite{end+1,1} = 30000;
bytesHeader = bytesHeader + 4;
precision{end+1,1} = 'bit32';
% 8. Time origin, 16 bytes, UTC time at which the data file is collected. 8
%       2-byte values definining Year, Month, DayOfWeek, Day, Hour, Minute,
%       Second, Millisecond
c = clock;
c = [c(1:2),1,c(3:end),0];
str = '';
for i = 1:numel(c)
    headerInfoWrite{end+1,1} = c(i);
    precision{end+1,1} = 'bit16';
end
bytesHeader = bytesHeader + 16;

% 9. Application to create file, 32 bytes, which program created the file,
%       null terminated
headerInfoWrite{end+1,1} = ['writeNev.m',char(0)];
while(length(headerInfoWrite{end}) <= 32)
    headerInfoWrite{end}(end+1) = ' ';
end
bytesHeader = bytesHeader + 32;
precision{end+1,1} = 'char';
% 10. Comment field, 256 bytes, null terminated string for comments
lenComments = length(comments);
if(lenComments > 255)
    comments = comments(1:255);
end
headerInfoWrite{end+1,1} = [comments,char(0)];

while(length(headerInfoWrite{end}) <= 254)
    headerInfoWrite{end,1}(end+1) = ' ';
end

bytesHeader = bytesHeader + 256;
precision{end+1,1} = 'char';
% 11. Number of extended headers, 4 bytes, a value indicating the number of
%       extended headers
headerInfoWrite{end+1,1} = numExtendedHeaders;
extendedHeaderIdx = numel(headerInfoWrite);
bytesHeader = bytesHeader + 4;
precision{end+1,1} = 'bit32';


%% Section 2 - header extend information
% this is used to add additional 32 byte fixed length additional
% configuration information (8 byte id, 24 byte information)
% 1. "NEUEVWAV" followed by:
%   a. ElectrodeID, 2 bytes, electrode ID number
%   b. Physical connector, 1 byte, physical connector (Front-End Bank A, B, C, D are
%       1, 2, 3, 4)
%   c. Connector Pin, 1 byte, physical system connector pin or channel
%       connected to electrode (1-37 on bank A, B, C, D)
%   d. Digitization scaling factor, 2 bytes, nV per LSB step (250?)
%   e. Energy threshold, 2 bytes, 0 if none used (0)
%   f. High threshold, 2 bytes, amplitude high threshold used (uV) (0)
%   g. Lowthreshold, 2 bytes, amplitude low threshold used (uV) (some #)
%   h. Number of sorted units, 1 byte, 0 for no unit classifcation
%   i. Bytes per waveform, 1 byte, number of bytes per waveform sample (2)
%   j. Spike Width (samples), 2 bytes, number of samples for each waveform
%       (default is 48 samples)
% 2. Others -- won't deal with until I have to later

% NEUEVWAV data stored in a mapfile -- read in mapfile
arrayMap=loadMapFile(mapfileName);
numExtendedHeaders = numExtendedHeaders + size(arrayMap,1);
for arrayMapIdx = 1:size(arrayMap,1)
    % write "NEUEVWAV"
    headerInfoWrite{end+1,1} = 'NEUEVWAV';
    precision{end+1,1} = 'char';
    % electrode ID -- 2 bytes
    headerInfoWrite{end+1,1} = arrayMap.chan(arrayMapIdx);
    precision{end+1,1} = 'bit16';
    % physical connector -- 1 byte (A,B,C,D = 1,2,3,4)
    headerInfoWrite{end+1,1} = arrayMap.bank{arrayMapIdx}-64;
    precision{end+1,1} = 'bit8';
    % connector pin -- 1 byte (1-37)
    headerInfoWrite{end+1,1} = arrayMap.pin(arrayMapIdx);
    precision{end+1,1} = 'bit8';
    % digitization scaling factor -- 2 bytes 
    headerInfoWrite{end+1,1} = 250;
    precision{end+1,1} = 'bit16';
    % Energy threshold, 2 bytes, 0 if none used (0)
    headerInfoWrite{end+1,1} = 0;
    precision{end+1,1} = 'bit16';
    % High threshold, 2 bytes, amplitude high threshold used (uV) (0)
    headerInfoWrite{end+1,1} = 0;
    precision{end+1,1} = 'bit16';
    % Lowthreshold, 2 bytes, amplitude low threshold used (uV) (some #)
    headerInfoWrite{end+1,1} = 0;
    precision{end+1,1} = 'bit16';
    % Number of sorted units, 1 byte, 0 for no unit classifcation
    headerInfoWrite{end+1,1} = 0;
    precision{end+1,1} = 'bit8';
    % Bytes per waveform, 1 byte, number of bytes per waveform sample (2)
    headerInfoWrite{end+1,1} = 2;
    precision{end+1,1} = 'bit8';
    % Spike Width (samples), 2 bytes, number of samples for each waveform
    headerInfoWrite{end+1,1} = 48;
    precision{end+1,1} = 'bit16';
    % reserved bytes, 8 bytes, write as 0
    headerInfoWrite{end+1,1} = 0;
    precision{end+1,1} = 'bit64';
end
headerInfoWrite{extendedHeaderIdx} = numExtendedHeaders;
headerInfoWrite{4,1} = bytesHeader; % this will be updated last

% NEUEVLBL
% packetID, 8 bytes, "NEUEVLBL"
% electrodeID, 2 bytes, electrode id (#)
% label, 16 bytes, label, null terminated
% remaining, 6 bytes, write as 0
numExtendedHeaders = numExtendedHeaders + size(arrayMap,1);
lbl2(1:8) = setstr(0); 
for arrayMapIdx = 1:size(arrayMap,1)
    % write "NEUEVLBL"
    headerInfoWrite{end+1,1} = 'NEUEVLBL';
    precision{end+1,1} = 'char';
    % electrode ID -- 2 bytes
    headerInfoWrite{end+1,1} = arrayMap.chan(arrayMapIdx);
    precision{end+1,1} = 'bit16';
    % electrode label -- 16 bytes, null terminated
    lbl1 = arrayMap.label{arrayMapIdx};
    while(length(lbl1)<8)
        lbl1(length(lbl1)+1)=setstr(0);
    end
    headerInfoWrite{end+1,1} = lbl1;
    precision{end+1,1} = 'char';
    headerInfoWrite{end+1,1} = lbl2;
    precision{end+1,1} = 'char';
    % reserved bytes, 6 bytes, write as 0
    headerInfoWrite{end+1,1} = 0;
    precision{end+1,1} = 'ubit48';
end

% NEUEVFLT
% packetID, 8 bytes, "NEUEVFLT"
% electrodeID, 2 bytes, electrode id (#)
% High freq corner, 4 bytes, in mHz
% High freq order, 4 bytes, 0 = NONE
% High filter type, 2 bytes, 0=none, 1=butter
% Low freq corner, 4 bytes, mHz
% Low freq order, 4 bytes, 0 = none
% Low filter type, 2 bytes, 0=none, 1=butter
% remaining, 2 bytes, write as 0

numExtendedHeaders = numExtendedHeaders + size(arrayMap,1);
for arrayMapIdx = 1:size(arrayMap,1)
    % write "NEUEVFLT"
    headerInfoWrite{end+1,1} = 'NEUEVFLT';
    precision{end+1,1} = 'char';
    % electrode ID -- 2 bytes
    headerInfoWrite{end+1,1} = arrayMap.chan(arrayMapIdx);
    precision{end+1,1} = 'bit16';
    % High freq corner, 4 bytes, in mHz
    headerInfoWrite{end+1,1} = 0;
    precision{end+1,1} = 'bit32';
    % High freq order, 4 bytes, 0 = NONE
    headerInfoWrite{end+1,1} = 0;
    precision{end+1,1} = 'bit32';
    % High filter type, 2 bytes, 0=none, 1=butter
    headerInfoWrite{end+1,1} = 0;
    precision{end+1,1} = 'bit16';
    % Low freq corner, 4 bytes, in mHz
    headerInfoWrite{end+1,1} = 0;
    precision{end+1,1} = 'bit32';
    % Low freq order, 4 bytes, 0 = NONE
    headerInfoWrite{end+1,1} = 0;
    precision{end+1,1} = 'bit32';
    % Low filter type, 2 bytes, 0=none, 1=butter
    headerInfoWrite{end+1,1} = 0;
    precision{end+1,1} = 'bit16';
    % reserved bytes, 2 bytes, write as 0
    headerInfoWrite{end+1,1} = 0;
    precision{end+1,1} = 'bit16';
end

% DIGLABEL
% packetID, 8 bytes, "DIGLABEL"
% label, 16 bytes, label, null terminated
% mode, 1 byte, 0 = serial, 1 = parallel
% remaining, 6 bytes, write as 0

numExtendedHeaders = numExtendedHeaders + 1;
% write "NEUEVFLT"
headerInfoWrite{end+1,1} = 'DIGLABEL';
precision{end+1,1} = 'char';
% label -- 16 bytes, null terminated
lbl1 = 'dummy';
while(length(lbl1)<8)
    lbl1(length(lbl1)+1)=setstr(0);
end
headerInfoWrite{end+1,1} = lbl1;
precision{end+1,1} = 'char';
headerInfoWrite{end+1,1} = lbl2;
precision{end+1,1} = 'char';
% mode, 1 byte, 0 = serial, 1 = parallel
headerInfoWrite{end+1,1} = 0;
precision{end+1,1} = 'bit8';
% reserved bytes, 7 bytes, write as 0
headerInfoWrite{end+1,1} = 0;
precision{end+1,1} = 'bit56';

bytesHeader = bytesHeader + numExtendedHeaders*32;
headerInfoWrite{extendedHeaderIdx,1} = numExtendedHeaders;
headerInfoWrite{4,1} = bytesHeader; % this will be updated last


%% write all of the header related information to the file

for idx = 1:numel(headerInfoWrite)
    fwrite(fid,headerInfoWrite{idx,1},precision{idx,1},0,'ieee-le');
end

%% Section 3 - data packets
% spike data
% each packet begins with a 4 byte time stamp and a 2 byte packed id
% 1. Packet IDs 1 through 2048 -- spike event, ID refers to elec number
%   a. Time stamp, 4 bytes, zero is beginning of data acquisition
%   b. Packet ID, 2 bytes, electrode number (1-96)
%   c. Unit Classification Number, 1 byte, 0 = unclassified
%   d. reserved, 1 byte, use 0
%   e. waveform, packetWidth-8 bytes, the sampled waveform of the spike
packetInfoWrite = {};
precision = {};
for spikeIdx = 1:numel(timeSpikes)
    % timestamp, 4 bytes, 0 is beginning of file
    fwrite(fid,floor(30000*timeSpikes(spikeIdx)),'bit32',0,'ieee-le');
%     packetInfoWrite{end+1,1} = floor(30000*timeSpikes(spikeIdx));
%     precision{end+1,1} = 'bit32';
    % packet id, 2 bytes, elec number
    fwrite(fid,chanSpikes(spikeIdx,1),'bit16',0,'ieee-le');
%     packetInfoWrite{end+1,1} = chanSpikes(spikeIdx,1);
%     precision{end+1,1} = 'bit16';
    % unit classification #, 1 byte, 0 is unclassified
    fwrite(fid,0,'bit16',0,'ieee-le');
%     packetInfoWrite{end+1,1} = 0;
%     precision{end+1,1} = 'bit8';
%     % reserved, 1 byte, 0
%     packetInfoWrite{end+1,1} = 0;
%     precision{end+1,1} = 'bit8';
    % waveform, packetWidth-8 bytes, sample waveform, 2 bytes per datapoint
    for waveformIdx = 1:size(waveformSpikes,2)
        fwrite(fid,waveformSpikes(spikeIdx,waveformIdx),'bit16',0,'ieee-le');
%         packetInfoWrite{end+1,1} = waveformSpikes(spikeIdx,waveformIdx);
%         precision{end+1,1} = 'bit16';
    end
end

% for idx = 1:numel(packetInfoWrite)
%     fwrite(fid,packetInfoWrite{idx,1},precision{idx,1},0,'ieee-le');
% end

% digital/serial inputs
% each packet begins with a 4 byte time stamp and a 2 byte packed id
% 1. Packet IDs 1 through 2048 -- spike event, ID refers to elec number
%   a. Time stamp, 4 bytes, zero is beginning of data acquisition
%   b. Packet ID, 2 bytes, electrode number (1-96)
%   c. Unit Classification Number, 1 byte, 0 = unclassified
%   d. reserved, 1 byte, use 0
%   e. waveform, packetWidth-8 bytes, the sampled waveform of the spike
packetInfoWrite = {};
precision = {};
digDummyTS = [0,timeSpikes(end)+1];
for dummyIdx = 1:numel(digDummyTS)
    % timestamp, 4 bytes, 0 is beginning of file
    packetInfoWrite{end+1,1} = floor(30000*digDummyTS(dummyIdx));
    precision{end+1,1} = 'bit32';
    % packet id, 2 bytes, 0
    packetInfoWrite{end+1,1} = 0;
    precision{end+1,1} = 'bit16';
    % packet insertion reason -- if digital channel then 1
    packetInfoWrite{end+1,1} = 1;
    precision{end+1,1} = 'bit8';
    % reserved, 1 byte, 0
    packetInfoWrite{end+1,1} = 0;
    precision{end+1,1} = 'bit8';
    % dig input, 2 bytes, 
    packetInfoWrite{end+1,1} = dummyIdx;
    precision{end+1,1} = 'bit16';
    % reserved, 10 bytes
    packetInfoWrite{end+1,1} = 0;
    precision{end+1,1} = 'bit64';
    packetInfoWrite{end+1,1} = 0;
    precision{end+1,1} = 'bit16';
end

for idx = 1:numel(packetInfoWrite)
    fwrite(fid,packetInfoWrite{idx,1},precision{idx,1},0,'ieee-le');
end


% close the NEV file
fclose(fid);



%% we also used to need to write a single NSx file to put it back into a cds later
% open the NSx file -- add extension if necessary
[pathstr,fname,fext]=fileparts(filename);
fid = fopen(strcat(fname,'.ns1'),'wb');

numExtendedHeaders = 0;
bytesHeader = 0;
headerInfoWrite = {};
precision = {};

%% Section 1 -- Basic Header
% 1. File type ID, 8 bytes, always set as "NEURALCD" for neural events
headerInfoWrite{1,1} = 'NEURALCD';
bytesHeader = bytesHeader + 8;
precision{1,1} = 'char';
% 2. File spec, 2 bytes, major and minor revision numbers -- 0x0203 for 2.3
headerInfoWrite{end+1,1} = 770;
bytesHeader = bytesHeader + 2;
precision{end+1,1} = 'int16';
% 3. Bytes in headers, 4 bytes -- total number of bytes in both headers.
%       Used as a zero idx reference to the first data packet
headerInfoWrite{end+1,1} = bytesHeader; % this will updated last
bytesHeader = bytesHeader + 4;
precision{end+1,1} = 'int32';
% 4. Label, 16 bytes, label of the sampling group, must be null terminated
temp(1:16)=setstr(0);
temp(1:7) = 'LFP Low';
headerInfoWrite{end+1,1} = temp; 
precision{end+1,1} = 'char';
% temp(1:8)=setstr(0);
% headerInfoWrite{end+1,1} = temp; 
% precision{end+1,1} = 'char';
bytesHeader = bytesHeader + 16;
% 5. Comment, 256 bytes, must be null terminated
temp(1:256)=setstr(0);
headerInfoWrite{end+1,1} = temp;
precision{end+1,1} = 'char';
bytesHeader = bytesHeader + 256;
% 6. Period, 4 bytes, number of 1/30000 between data points
headerInfoWrite{end+1,1} = 60;
precision{end+1,1} = 'int32';
bytesHeader = bytesHeader + 4;
% 7. Time resolution of samples, 4 bytes, sampling frequency used to
%       digitize neural waveforms. Will be set to 30,000
headerInfoWrite{end+1,1} = 30000;
bytesHeader = bytesHeader + 4;
precision{end+1,1} = 'int32';
% 8. Time origin, 16 bytes, UTC time at which the data file is collected. 8
%       2-byte values definining Year, Month, DayOfWeek, Day, Hour, Minute,
%       Second, Millisecond
c = clock;
c = [c(1:2),1,c(3:end),0];
for i = 1:numel(c)
    headerInfoWrite{end+1,1} = c(i);
    precision{end+1,1} = 'int16';
end
bytesHeader = bytesHeader + 16;

% 9. Channel Count, 4 bytes, number of channels per data point
headerInfoWrite{end+1,1} = 1;
bytesHeader = bytesHeader + 4;
precision{end+1,1} = 'int32';

%% Extended header 1
numExtendedHeaders = numExtendedHeaders + 1;
% write "CC", 2 bytes
% headerInfoWrite{end+1,1} = 67; % 'C'
% precision{end+1,1} = 'int8';
% headerInfoWrite{end+1,1} = 67; % 'C'
% precision{end+1,1} = 'int8';
headerInfoWrite{end+1,1} = 'CC'; % 'C'
precision{end+1,1} = 'char';
% electrode ID -- 2 bytes
headerInfoWrite{end+1,1} = 97;
precision{end+1,1} = 'int16';
% label, 16 bytes, label of the electrode, must be null terminated
temp = '';
temp(1:16)=setstr(0);
temp(1:7) = 'elecJOE';
headerInfoWrite{end+1,1} = temp; 
precision{end+1,1} = 'char';
% physical connector -- 1 byte (A,B,C,D = 1,2,3,4)
headerInfoWrite{end+1,1} = 1;
precision{end+1,1} = 'uint8';
% connector pin -- 1 byte (1-37)
headerInfoWrite{end+1,1} = 1;
precision{end+1,1} = 'uint8';
% min digital value
headerInfoWrite{end+1,1} = -1000;
precision{end+1,1} = 'int16';
% max digital value
headerInfoWrite{end+1,1} = 1000;
precision{end+1,1} = 'int16';
% min analog value
headerInfoWrite{end+1,1} = -1000;
precision{end+1,1} = 'int16';
% max analog value
headerInfoWrite{end+1,1} = 1000;
precision{end+1,1} = 'int16';
% units, 16 bytes, label of the electrode, must be null terminated
temp(1:16) = setstr(0);
temp(1:3)= 'JOE';
headerInfoWrite{end+1,1} = temp; % this will updated last
precision{end+1,1} = 'char';
% High freq corner, 4 bytes, in mHz
headerInfoWrite{end+1,1} = 0;
precision{end+1,1} = 'int32';
% High freq order, 4 bytes, 0 = NONE
headerInfoWrite{end+1,1} = 0;
precision{end+1,1} = 'int32';
% High filter type, 2 bytes, 0=none, 1=butter
headerInfoWrite{end+1,1} = 0;
precision{end+1,1} = 'int16';
% Low freq corner, 4 bytes, in mHz
headerInfoWrite{end+1,1} = 0;
precision{end+1,1} = 'int32';
% Low freq order, 4 bytes, 0 = NONE
headerInfoWrite{end+1,1} = 0;
precision{end+1,1} = 'int32';
% Low filter type, 2 bytes, 0=none, 1=butter
headerInfoWrite{end+1,1} = 0;
precision{end+1,1} = 'int16';

bytesHeader = bytesHeader + numExtendedHeaders*66; % this is updated below
headerInfoWrite{3,1} = bytesHeader; 
%% write all of the header related information to the file

for idx = 1:numel(headerInfoWrite)
    fwrite(fid,headerInfoWrite{idx,1},precision{idx,1},0,'ieee-le');
end

%% Section 3 - data packets
packetInfoWrite = {};
precision = {};
for idx = 1:2
    % header, 1 byte, 0x01
    packetInfoWrite{end+1,1} = 1;
    precision{end+1,1} = 'int8';
    % timestamp, 4 bytes, 0 is beginning of file
    if(idx==1)
        packetInfoWrite{end+1,1} = 0.1;
    else
        packetInfoWrite{end+1,1} = ceil(30000*timeSpikes(end))+1;
    end
    precision{end+1,1} = 'int32';
    % # data points, 4 bytes,
    packetInfoWrite{end+1,1} = 1;
    precision{end+1,1} = 'int32';
    % data point, 2 bytes
    packetInfoWrite{end+1,1} = idx*10;
    precision{end+1,1} = 'int16';
end

for idx = 1:numel(packetInfoWrite)
    fwrite(fid,packetInfoWrite{idx,1},precision{idx,1},0,'ieee-le');
end

% close the NSx file
fclose(fid);


end




% %% we also used to need to write a single NSx file to put it back into a cds later
% % open the NSx file -- add extension if necessary
% [pathstr,fname,fext]=fileparts(filename);
% fid = fopen(strcat(fname,'.ns1'),'wb');
% 
% numExtendedHeaders = 0;
% bytesHeader = 0;
% headerInfoWrite = {};
% precision = {};
% 
% %% Section 1 -- Basic Header
% % 1. File type ID, 8 bytes, always set as "NEURALCD" for neural events
% headerInfoWrite{1,1} = 'NEURALCD';
% bytesHeader = bytesHeader + 8;
% precision{1,1} = 'char';
% % 2. File spec, 2 bytes, major and minor revision numbers -- 0x0203 for 2.3
% headerInfoWrite{end+1,1} = 770;
% bytesHeader = bytesHeader + 2;
% precision{end+1,1} = 'int16';
% % 3. Bytes in headers, 4 bytes -- total number of bytes in both headers.
% %       Used as a zero idx reference to the first data packet
% headerInfoWrite{end+1,1} = bytesHeader; % this will updated last
% bytesHeader = bytesHeader + 4;
% precision{end+1,1} = 'int32';
% % 4. Label, 16 bytes, label of the sampling group, must be null terminated
% temp(1:16)=setstr(0);
% temp(1:7) = 'LFP Low';
% headerInfoWrite{end+1,1} = temp; 
% precision{end+1,1} = 'char';
% % temp(1:8)=setstr(0);
% % headerInfoWrite{end+1,1} = temp; 
% % precision{end+1,1} = 'char';
% bytesHeader = bytesHeader + 16;
% % 5. Comment, 256 bytes, must be null terminated
% temp(1:256)=setstr(0);
% headerInfoWrite{end+1,1} = temp;
% precision{end+1,1} = 'char';
% bytesHeader = bytesHeader + 256;
% % 6. Period, 4 bytes, number of 1/30000 between data points
% headerInfoWrite{end+1,1} = 60;
% precision{end+1,1} = 'int32';
% bytesHeader = bytesHeader + 4;
% % 7. Time resolution of samples, 4 bytes, sampling frequency used to
% %       digitize neural waveforms. Will be set to 30,000
% headerInfoWrite{end+1,1} = 30000;
% bytesHeader = bytesHeader + 4;
% precision{end+1,1} = 'int32';
% % 8. Time origin, 16 bytes, UTC time at which the data file is collected. 8
% %       2-byte values definining Year, Month, DayOfWeek, Day, Hour, Minute,
% %       Second, Millisecond
% c = clock;
% c = [c(1:2),1,c(3:end),0];
% for i = 1:numel(c)
%     headerInfoWrite{end+1,1} = c(i);
%     precision{end+1,1} = 'int16';
% end
% bytesHeader = bytesHeader + 16;
% 
% % 9. Channel Count, 4 bytes, number of channels per data point
% headerInfoWrite{end+1,1} = 1;
% bytesHeader = bytesHeader + 4;
% precision{end+1,1} = 'int32';
% 
% %% Extended header 1
% numExtendedHeaders = numExtendedHeaders + 1;
% % write "CC", 2 bytes
% % headerInfoWrite{end+1,1} = 67; % 'C'
% % precision{end+1,1} = 'int8';
% % headerInfoWrite{end+1,1} = 67; % 'C'
% % precision{end+1,1} = 'int8';
% headerInfoWrite{end+1,1} = 'CC'; % 'C'
% precision{end+1,1} = 'char';
% % electrode ID -- 2 bytes
% headerInfoWrite{end+1,1} = 97;
% precision{end+1,1} = 'int16';
% % label, 16 bytes, label of the electrode, must be null terminated
% temp = '';
% temp(1:16)=setstr(0);
% temp(1:7) = 'elecJOE';
% headerInfoWrite{end+1,1} = temp; 
% precision{end+1,1} = 'char';
% % physical connector -- 1 byte (A,B,C,D = 1,2,3,4)
% headerInfoWrite{end+1,1} = 1;
% precision{end+1,1} = 'uint8';
% % connector pin -- 1 byte (1-37)
% headerInfoWrite{end+1,1} = 1;
% precision{end+1,1} = 'uint8';
% % min digital value
% headerInfoWrite{end+1,1} = -1000;
% precision{end+1,1} = 'int16';
% % max digital value
% headerInfoWrite{end+1,1} = 1000;
% precision{end+1,1} = 'int16';
% % min analog value
% headerInfoWrite{end+1,1} = -1000;
% precision{end+1,1} = 'int16';
% % max analog value
% headerInfoWrite{end+1,1} = 1000;
% precision{end+1,1} = 'int16';
% % units, 16 bytes, label of the electrode, must be null terminated
% temp(1:16) = setstr(0);
% temp(1:3)= 'JOE';
% headerInfoWrite{end+1,1} = temp; % this will updated last
% precision{end+1,1} = 'char';
% % High freq corner, 4 bytes, in mHz
% headerInfoWrite{end+1,1} = 0;
% precision{end+1,1} = 'int32';
% % High freq order, 4 bytes, 0 = NONE
% headerInfoWrite{end+1,1} = 0;
% precision{end+1,1} = 'int32';
% % High filter type, 2 bytes, 0=none, 1=butter
% headerInfoWrite{end+1,1} = 0;
% precision{end+1,1} = 'int16';
% % Low freq corner, 4 bytes, in mHz
% headerInfoWrite{end+1,1} = 0;
% precision{end+1,1} = 'int32';
% % Low freq order, 4 bytes, 0 = NONE
% headerInfoWrite{end+1,1} = 0;
% precision{end+1,1} = 'int32';
% % Low filter type, 2 bytes, 0=none, 1=butter
% headerInfoWrite{end+1,1} = 0;
% precision{end+1,1} = 'int16';
% 
% bytesHeader = bytesHeader + numExtendedHeaders*66; % this is updated below
% headerInfoWrite{3,1} = bytesHeader; 
% %% write all of the header related information to the file
% 
% for idx = 1:numel(headerInfoWrite)
%     fwrite(fid,headerInfoWrite{idx,1},precision{idx,1},0,'ieee-le');
% end
% 
% %% Section 3 - data packets
% packetInfoWrite = {};
% precision = {};
% for idx = 1:2
%     % header, 1 byte, 0x01
%     packetInfoWrite{end+1,1} = 1;
%     precision{end+1,1} = 'int8';
%     % timestamp, 4 bytes, 0 is beginning of file
%     if(idx==1)
%         packetInfoWrite{end+1,1} = 0.1;
%     else
%         packetInfoWrite{end+1,1} = ceil(30000*timeSpikes(end));
%     end
%     precision{end+1,1} = 'int32';
%     % # data points, 4 bytes,
%     packetInfoWrite{end+1,1} = 1;
%     precision{end+1,1} = 'int32';
%     % data point, 2 bytes
%     packetInfoWrite{end+1,1} = idx*10;
%     precision{end+1,1} = 'int16';
% end
% 
% for idx = 1:numel(packetInfoWrite)
%     fwrite(fid,packetInfoWrite{idx,1},precision{idx,1},0,'ieee-le');
% end
% 
% % close the NSx file
% fclose(fid);

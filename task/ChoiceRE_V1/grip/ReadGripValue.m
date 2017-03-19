function [Values, Times] = ReadGripValue(DeviceName, Handle, Channels)
% Read value on grip device.
%
% Usage: [Values, Times] = ReadGripValue(DeviceName, Handle, Channels)
%
% DeviceName: 'UsbCED' for the CED (old) device
%             'SerialMBB' for the (new) device using the serial port
%             'MIE' for the pinch grip my MIE medical research ltd
% Channels: [1], [2] or [1 2] 
%           With 1 (/2): the channel label 1 (/2) on the SerialMBB device 
%           and 0 (/1) on the CED ADC input
% Handle: the device handle
%
% Values: the returned value (if several Channels, they are returned in the 
%         same order as requested)
% Times: Time point (PsychToolBox GetSecs.m) for each measurement.
%
%   See also: InitializeGrip, ReadGripValue
%
% Florent Meyniel 2013-02-13
% 

Values = zeros(size(Channels));
Times  = zeros(size(Channels));
    
if strcmpi(DeviceName, 'UsbCED')

    Channels = Channels - 1 ; % coded as 0 and 1

    matced32('cedSendString','adcmem,?;'); % ask the 1401 about its state
    while(matced32('cedStat1401')==0);
        drawnow; % EF 9-11-04
    end;

    wh = str2num(matced32('cedGetString'));	% which half is filling?
    if (wh == 2)
        r = matced32('cedToHost',500,0)';
    else
        r = matced32('cedToHost',500,1000)';
    end

    for iChan = 1:length(Channels)
        Values(iChan) = nanmean(r(iChan:length(Channels):end));
        Times(iChan) = GetSecs;
    end

elseif strcmpi(DeviceName, 'SerialMBB')

    for iChan = 1:numel(Channels)
        if     Channels(iChan) == 1
            IOPort('Write', Handle, uint8('A'));
        elseif Channels(iChan) == 2
            IOPort('Write', Handle, uint8('C'));
        else    error('Unknown Channel %d, should be 1 or 2', Channels)
        end
        WaitSecs(0.005); % wait for Microcontroler's answer
        [data, Times(iChan)] = IOPort('Read', Handle);
        Values(iChan) = str2double(char(data));
    end
    
elseif strcmpi(DeviceName, 'MIE')
    
    tmp = [];
    IOPort('Purge', Handle);
    t = GetSecs;
    while length(tmp)~=2 && (GetSecs - t) < 2
        WaitSecs(0.02);
        [tmp, Times] = IOPort('Read', Handle);
    end
    
    if GetSecs - t > 2
        Values = NaN;
        Times = NaN;
    else
        Byte1 = dec2bin(tmp(1), 8);
        Byte2 = dec2bin(tmp(2), 8);
        
        % INFORMATION BY THE CONSTRUCTOR:
        % Buffer comprises of two bytes:
        %
        %				 BYTE 1	  BYTE 2
        %	BIT			87654321 87654321
        %	PURPOSE		DDDDxxSx DDDDDDDD (where D=DATA, S=SIGN, x=DON'T CARE)
        %  REPRESENTS	CBA9  S	 43218765 (bit number of value in hex - base 16)
        %
        %	Byte 1 : Bit 2 is the sign bit (indicated by 'S')
        %	The 12 bit value is calculated as follows:
        %		Multiply 4 most significant bits of byte 1 by 256.
        %		Add 4 most significant bits of byte 2 shifted right 4.
        %		Add 4 least significant bits of byte 2 shifted left 4.
        %
        %	Thus Value = (((Byte1 & 240) >> 4) * 256) |
        %	  ((Byte2 & 240) >> 4) | ((Byte2 & 15) << 4);
        %
        %	Then deal with negation:
        %		if(!(Byte1 & 2))		if bit 2 of byte 1
        %			Value *= -1;
        %
        % Convert to volts using channel range (default +5V to -5V = 5 - -5 = 10)
        % Assumes 12 bits (8191)
        % Subtract zero offset measured before test commenced
        % double dVolts = ((double)iCasMap * (double)m_fChannelRange / (double)8191) - m_fChannelZeroOffset;
        %
        % Convert to newtons
        % Uses calibration value for digital analyser (default for standard units: 333.0 = 0.333 Volts)
        % For myometer uses joint distance to strap in mm / 1000 (result in Newton Metres)
        % For pinchgrip uses a joint distance of 1.0 (result in Newtons)
        % double dNewtons = (double)(dVolts * m_fChannelCalibrationValue * m_fJointDistance);
        
        
        if Byte1(7); Values = bin2dec([Byte1(1:4) Byte2(5:8) Byte2(1:4)]);
        else Values = -1 * bin2dec([Byte1(1:4) Byte2(5:8) Byte2(1:4)]); end
        
        % convert to newtons (according to the constructor)
        Values = Values * 10 / 8191 * 333;
    end

else
    error('Unknow method %s. It should be ''UsbCED'', ''MIE'' or ''SerialMBB''', DeviceName)
end


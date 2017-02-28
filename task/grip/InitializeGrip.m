function Handle = InitializeGrip(DeviceName, Channels)
% Initialize the grip device
%
% Usage: Handle = InitializeGrip(DeviceName, Channels)
%
% DeviceName: 'UsbCED' for the CED (old) device
%             'SerialMBB' for the (new) device using the serial port
% Channels: [1], [2] or [1 2] 
%           With 1 (/2): the channel label 1 (/2) on the SerialMBB device
%           and 0 (/1) on the CED ADC input
%
% Caution: depending on the CED and / or Matlab version, some path issue
% may arise... see in the script.
%
%   See also: ReadGripValue, CloseGripDevice
%
% Florent Meyniel 2013-02-13
% 

Handle = [];

% Check that Psychtoolbox is installed
try PsychtoolboxVersion
catch
    error('PSYCHTOOLBOX IS NOT INSTALLED!')
end

if strcmpi(DeviceName, 'UsbCED')
    
    Channels = Channels - 1 ; % coded as 0 and 1

    % Open CED
    matced32('cedOpen');
    
    % Solution 1: 1401Lang is not in the path: it should be specified when
    % calling the DLL:
%     LANG1401PATH = 'C:\1401Lang\';
%     if ~exist(LANG1401PATH, 'dir'); error('Cannot find Lang1401 at %s', LANG1401PATH); end
%     matced32('cedLdX', LANG1401PATH, 'ADCMEM'); % 2nd argument: the full path
%     
    % Soluation 2: if 1401Lang is in Matlab's paths: 
    matced32('cedLd', 'ADCMEM');
    
    % Specify how to read data
    matced32('cedLd', 'event');
    
    % Initialise the settings
    matced32('cedSendString', 'event, d, 26;');
    matced32('cedSendString', 'event, m, 128;');
    
    % Sets up transfer from ADC to memory
    % - i: interrupt driven,
    % 2: 2 bytes per point,
    % 0: mem start byte,
    % 2000: bytes in buffer,
    % {0 1 2 3}: inputs for each gripper,
    % 0: rpt,
    % ct: wait for ttl on E4,
    % 2: prescale,
    % 10: count 
    matced32( 'cedSendString', ['adcmem, i, 2, 0, 2000, ' num2str(Channels) ', 0, ct, 2, 10;'] );
    matced32('cedSendString', 'event, i, 26;');

elseif strcmpi(DeviceName, 'SerialMBB')
    
    % Serial Port Name Usually: 'COM1' on Windows; '/dev/ttyS0' on Linux
    % Change if necessary
    if ispc
        PortName = 'COM1';
    end
    if isunix
        PortName = '/dev/ttyS0';
    end
    if ismac
        error('Mac User! should specify the Serial port address in the script')
    end

    % intialize the serial port
    configString = strcat('BaudRate=38400', ' ', ...
        'Parity=None', ' ', ...
        'DataBits=8', ' ', ...
        'StopBits=1');

    Handle = IOPort('OpenSerialPort', PortName, configString);
    IOPort('Purge', Handle);
    
    % read a value, to be discarded (as on some computers, the first value
    % read is NaN...)
    for iChan = 1:numel(Channels)
        if     Channels(iChan) == 1
            IOPort('Write', Handle, uint8('A'));
        elseif Channels(iChan) == 2
            IOPort('Write', Handle, uint8('C'));
        else    error('Unknown Channel %d, should be 1 or 2', Channels)
        end
        WaitSecs(0.005); % wait for Microcontroler answer
        IOPort('Read', Handle);
    end
    
elseif strcmpi(DeviceName, 'MIE')
    
    % Serial Port Name Usually: 'COM1' on Windows; '/dev/ttyS0' on Linux
    % Change if necessary
    if ispc
        PortName = 'COM1';
    end
    if isunix
        PortName = '/dev/ttyS0';
    end
    if ismac
        error('Mac User! should specify the Serial port address in the script')
    end
    
    % intialize the serial port
    configString = strcat('BaudRate=9600', ' ', ...
    'Parity=None', ' ', ...
    'DataBits=8', ' ', ...
    'StopBits=2', ' ', ...
    'FlowControl=None');

    Handle = IOPort('OpenSerialPort', PortName, configString);
    IOPort('Purge', Handle);    
    
else
    error('Unknow method %s. It should be ''UsbCED'', ''MIE'' or ''SerialMBB''', DeviceName)
end


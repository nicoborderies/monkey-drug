function CloseGripDevice(DeviceName, Handle)
% Close the grip device.
%
% Usage: CloseGripDevice(DeviceName, Handle)
%
% See also InitializeGrip, ReadGripValue
%
% Florent Meyniel 2013-02-13
% 


if strcmpi(DeviceName, 'UsbCED')

    matced32('cedReset');
    clear matced32.dll;
    clear use1432.dll;
    matced32('cedClose');

elseif strcmpi(DeviceName, 'SerialMBB') || strcmpi(DeviceName, 'MIE')

    IOPort('Close', Handle)
    
else
    error('Unknow method %s. It should be ''UsbCED'', ''MIE'' or ''SerialMBB''', DeviceName)
end


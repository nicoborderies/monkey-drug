function config_io

global cogent;

%create IO64 interface object
clear io64;
cogent.io.ioObj = io64;

%install the inpout32.dll driver
%status = 0 if installation successful
cogent.io.status = io64(cogent.io.ioObj);
if(cogent.io.status ~= 0)
    disp('inpout64 installation failed!')
else
    disp('inpout64 (re)installation successful.')
end


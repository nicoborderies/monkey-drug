%% test_fluid_distributor

% initialize
portdir = 'K:\monkeyResponse\IO32';
addpath(genpath(portdir));
config_io;
outp(888,0); % put all pins to zero

% parameters
Device.reward.time2volume = 0.9;
Device.reward.volumeUnit = 1.5;
Device.reward.timeUnit = Device.reward.volumeUnit/Device.reward.time2volume;
Device.reward.pauseDuration = 1.5;

% testing loop
clc;disp('Press any key to exit.');
exit=0;
KbQueueCreate;KbQueueStart;
while exit==0
    
    outp(888,00000001);
    WaitSecs(Device.reward.timeUnit);
    outp(888,0);
    WaitSecs(Device.reward.pauseDuration);
    keypress = KbQueueCheck;
    if keypress
        exit=1;
    end
end

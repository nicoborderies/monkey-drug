% Script to mimic an oscilloscope in Matlab for the grip force devices.
% This version includes the compatibility wrapper to works with both the
% 'UsbCED' or 'SerialMBB' grip devices.
% The script open a window were the data are displayed.
% The offset, scaling, sampling rate and channels can be modified by 
% calling a dialog box after pressing the 'c' key.
% To quit the script, press the escape key.
%
% Florent Meyniel 2012-01-03

% SET DEFAULTS VALUES
% ===================

% add grip device compatibility wrapper
p = mfilename('fullpath');
if isunix; ind = strfind(p, '/'); end
if ispc; ind = strfind(p, '\'); end
pathdir = p(1:ind(end));
% addpath(strcat(pathdir, 'GripCompatFunc'));

DeviceName = 'SerialMBB'; % 'UsbCED', 'SerialMBB' or 'MIE'
Channels = [1];
Handle = InitializeGrip(DeviceName, Channels);

wh_px    = 400;   % window's heigh
ww_px    = 400;   % window's width
VDiv     = 1000;  % vertical scaling ('voltage')
Offset   = 900;   % offset

% INITIALIZE DISPLAY
% ==================
KbName('UnifyKeyNames');
% get screen resolution
[w_px, h_px] = Screen('WindowSize', 0);
Fs = ceil(Screen('FrameRate', 0));

% check resolution
if h_px < wh_px
    error('Screen resolution mismatch (%d for a %d window)', h_px, wh_px)
end
if w_px < ww_px
    error('Screen resolution mismatch (%d for a %d window)', w_px, ww_px)
end

wd_osc = Screen('OpenWindow', 0, [], [1 1 1+wh_px 1+ww_px]);

fprintf(strcat('\n\n \t\t === MATOS ===\n', ...
               '\n   The channel, offset, scaling and sampling rate\n', ...
                 '   can be modified: call a dialog box \n', ...
                 '   by pressing the ''c'' key.\n\n'  ,...
                 '   To quit the script, press the escape key.\n\n'))

% LAUNCH OSCILLOSCOPE MODE
% ========================
% initialize variables
sampled = zeros(1,ww_px);
exit    = 0;
t       = GetSecs;
prevT   = 0;
mFs     = zeros(1, 10); % to estimate mean frequency sampling
while exit==0    
    % to escape the program
    [keyisdown, secs, keycode] = KbCheck;
    if keyisdown==1 && keycode(KbName('ESCAPE'))
        exit = 1;
    end
    
    % get voltage estimate
    [Values, Times] = ReadGripValue(DeviceName, Handle, Channels);
    sampled = [sampled(2:end), Values];
    
    eFs = 1/(Times(1) - prevT); % instantaneous true sampling frequency
    prevT = Times;
    mFs = [eFs mFs(1:9)];       % 10 last true sampling frequnecies 
    
    % get sampling rate, scale and offset
    if keyisdown==1 && keycode(KbName('c'))
       [Fs, VDiv, Offset, Channels] = getUserParam(Fs, VDiv, Offset, Channels);
    end
    
    % plot into a window
    xy = [1:ww_px; sampled];
    PlotDataOnOsc(wd_osc, xy, Offset, VDiv, wh_px, ww_px)
    
    % display parametersc
    DispParam(wd_osc, Offset, VDiv, xy(2,end), eFs, mean(mFs))
    
    % get timing, flush event and impose a sampling rate
    Screen('Flip', wd_osc);
    t = WaitSecs('UntilTime', t(end)+1/Fs);
end

% return to the original display mode
Screen('CloseAll');

CloseGripDevice(DeviceName, Handle)

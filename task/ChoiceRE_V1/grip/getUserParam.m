function [Fs, VDiv, Offset, channel] = getUserParam(Fs, VDiv, Offset, channel)
% get new parameters


prompt              = {'Channel' 'Offset' 'VDiv' 'Sampling rate'};
name                = 'Set New Parameters';
nLines              = 1;
defaultAnswer       = {num2str(channel) num2str(Offset) num2str(VDiv) num2str(Fs)};
options.WindowStyle = 'modal';
options.Interpreter = 'none';

answer = inputdlg(prompt, name, nLines, defaultAnswer, options);

if ~isempty(strmatch('This functionality is no longer supported under the -nojvm startup option', lastwarn))
    fprintf('\n Sorry... the dialog box is not supported with nojvm! type in the command line window\n')
    channel = input('Channel: ');
    Offset  = input('Offset: ');
    VDiv    = input('VDiv: ');
    Fs      = input('Fs: ');
else
   
    Offset  = str2double(answer{2});
    VDiv    = str2double(answer{3});
    Fs      = str2double(answer{4});
    channel = str2double(answer{1});
    
end





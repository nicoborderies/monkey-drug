%% test_dictaphone

% load
soundfile = 'xylophone_sound.wav';
[s,f] = audioread(soundfile);s = s';


% init
InitializePsychSound(1);
Device.audio.freq = 1000;
Device.audio.duration = 0.10;
Device.audio.intensity = 1;
Device.audio.handle = PsychPortAudio('Open', [], 1, 1, [], 2);
PsychPortAudio('Volume', Device.audio.handle, Device.audio.intensity);
Device.audio.beep = MakeBeep(Device.audio.freq, Device.audio.duration);

sound1 = PsychPortAudio('CreateBuffer', Device.audio.handle, [Device.audio.beep; Device.audio.beep]);
sound2 = PsychPortAudio('CreateBuffer', Device.audio.handle, s);

% play
PsychPortAudio('FillBuffer', Device.audio.handle, sound1);
PsychPortAudio('FillBuffer', Device.audio.handle, sound2);
PsychPortAudio('Start', Device.audio.handle, 1, 0, 0);

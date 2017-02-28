function DispParam(wd_osc, Offset, VDiv, val, eFs, mFs)
% display parameters

txt1 = sprintf(' Offset: %d', Offset);
txt2 = sprintf('   VDiv: %d', VDiv);
txt3 = sprintf('   curr: %5.0f', val);
txt4 = sprintf('true Fs: %5.0f', eFs);
txt5 = sprintf('mean Fs: %5.0f', mFs);

Screen('DrawText', wd_osc, txt1, 0, 0);
Screen('DrawText', wd_osc, txt2, 0, 20);
Screen('DrawText', wd_osc, txt3, 0, 40);
Screen('DrawText', wd_osc, txt4, 0, 60);
Screen('DrawText', wd_osc, txt5, 0, 80);


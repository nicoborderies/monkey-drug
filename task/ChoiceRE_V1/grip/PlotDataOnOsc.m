function PlotDataOnOsc(wd_osc, xy, Offset, VDiv, wh_px, ww_px)
% plot data on the oscilloscope window

% scale the display
SizeDiv = round(wh_px/5);
PlotXY(1,:) = xy(1,:);
PlotXY(2,:) = round(wh_px/2 - (xy(2,:)-Offset)*SizeDiv/VDiv);
Screen('DrawDots', wd_osc, PlotXY, 5, [0 0 0], [], 2);


% draw Oscillo reference lines
for i = 0:4 
Screen('DrawLine', wd_osc, [200 200 200],...
    0, ...
    round(wh_px/10) + i*SizeDiv, ...
    ww_px, ...
    round(wh_px/10) + i*SizeDiv);
end



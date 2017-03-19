

%% polynomial fit
% [p,s,mu] = polyfit(t,d,6);
% f_y = polyval(p,t,[],mu);

%% moving average
f_y = smooth(d,'moving'); % methods: 'loess' / 'sgolay' / 'rloess'


d2 = d - f_y;
f2 = f + d2;

figure; hold on
plot(d);
plot(f_y);
plot(d2);
plot(f);
plot(f2);
function [ z , w ] = smoothForceSignal( x , y )

 % time-series preprocessing: temporal smoothing / detect premature responses / extract temporal derivative
 
        x = x{:};
        y = y{:};
 
        window=2; 
        h=ones(window,1)/window;
            
            
        trialDuration = y(end)-y(1);
        dt=(trialDuration)/length(x);
        
        offset = x(1);
%         x = x-offset;
        
        Fdynamic_temp=filter(h,1,x);
        Fdynamic=fliplr(filter(h,1,fliplr(Fdynamic_temp)));
        z = Fdynamic;
        
        
        speed=diff(Fdynamic);
        Fspeed_temp=filter(h,1,speed);
        Fspeed=fliplr(filter(h,1,fliplr(Fspeed_temp)))/dt;
        w= Fspeed;
        
%     % plot signal smoothign
%     clf('reset');
%     hold on;
%     h1 = plot(y,x,'k');
%     h2 = plot(y,z,'Color',[0.5 0.5 0.5]);
%     h3 = plot(y(2:end),w,'Color','r');
%     xlabel('time');
%     ylabel('signal');
%     title('smoothing');
%     legend([h1,h2,h3],{'raw signal','filtered signal','derivative'})
%     hold off;
%     pause


end


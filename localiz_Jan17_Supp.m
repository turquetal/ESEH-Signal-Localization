clear all
close all
clc
load 2bar_acceldata_161213.mat
% load GV.mat

addpath C:\Users\usr1\Dropbox\Super_Mega_Script_Folder
addpath C:\Users\Dhareion\Dropbox\Super_Mega_Script_Folder
addpath C:\Users\User1\Dropbox\Super_Mega_Script_Folder
addpath /Users/turquetal/Dropbox/Super_Mega_Script_Folder
% % % % % %  Event Tracking over Injection Experimental Data 
evdetect=0;
arrivcheck=1;

if evdetect
thHF=4.3; % threshold high freq
taperzone=0; %sec 
dt=1e-6; 
f1=1e2; 
f2=1e4; % 1e4 for low fr events
% f2=4.5e5;
fHF1=1e5;
fHF2=4e5;
wSTA=5e-5;%sec
wLTA=1e-3; %sec
% wSTA=1e-4;%sec
% wLTA=1e-3; %sec
slide=wSTA; % jump of the sliding windows, sec 
fig=1; % if fig=1 it plots each file with high and low freq picks else it runs without plotting


for i=1:4
data(:,i)=taper(bp(data(:,i),dt,2e2,4.5e5),taperzone/dt);
end
nsta=wSTA/dt;
nlta=wLTA/dt;
%nw=floor((3600-wLTA-2*taperzone)/slide);
%t=wLTA/2+taperzone+[0:slide:(nw-1)*slide];
% for i=1:length(dat) % loop on the data  
    
%     a=char(dat(i));
%     %day and hour
%     day=(str2num(a(6:6))-6)*30+(str2num(a(7:8))-18);
%     hh=str2num(a(10:11));
    %read - filter data 
%     X0=bp(read_sac(char(dat(i))),0.01,0.01,40); %read data
% X=[data(inf:sup,1);postdata(:,1)];
X=[data(:,1)];
%     X=taper(rmean(X0),taperzone/dt);
    nw=floor((length(X)*dt-wLTA-2*taperzone)/slide);
    t=wLTA/2+taperzone+[0:slide:(nw-1)*slide];
%     char(dat(i))
    Y=taper(bp(X,dt,f1,f2),taperzone/dt); %filter data low freq
    YHF=taper(bp(X,dt,fHF1,fHF2),taperzone/dt);
    oth=0;
    oth_HF=0;
    oth_energy=0;
    
    T=zeros(size(t));
    THF=zeros(size(t));
    clear Z;
    cntr=1;
    cntrs=1;
    cntre=1;
    Z_Energy(1)=0;
    for j=1:nw % loop on the windows 
       t1=(slide*(j-1))/dt+floor((wLTA-wSTA)/(2*dt))+floor(taperzone/dt); %define t1 of small windows 
       t2=t1+floor(wSTA/dt);
%        
%        t2_1=1+floor((slide*(j-1))/dt+floor(taperzone/dt)); %only STA values
%        t2_2=t1+floor(wSTA/dt);
%        range=
%        Z_sta(j)=(sum(Y(t2_1:t2_2).^2)); %STA/LTA for low frequency 
%          [wavein0,f0]=fftrl(Y(t2_1:t2_2),time(t2_1:t2_2),0);
%      exist powini;
%          if ans
%              pownew=abs(wavein0).^2*dt;
%              Z_Energy(j)=trapz(abs(log(pownew)-log(powini)))/trapz(log(powini));
%                 powini=pownew;    
%          else
%              powini=abs(wavein0).^2*dt;
%          end
       
        tL1=floor((slide*(j-1))/dt+1+floor(taperzone/dt)); %define t1 of Large windows 
        tL2=tL1+floor(wLTA/dt)-1;
%         range=[[tL1:floor(t1-nsta/2)],[floor(t2+nsta/2):tL2]];
        range=[tL1:tL2];
             
%        Z(j)=sqrt(sum(Y(t1:t2)'*Y(t1:t2)))/sqrt(sum(Y(tL1:tL2)'*Y(tL1:tL2))); 
%                  ZHF(j)=(sum(YHF(t1:t2).*YHF(t1:t2)))/(sum(YHF(range).*YHF(range)))*length(range)/nsta; %STA/LTA for high frequency 
        ZHF(j)=(sum(YHF(t1:t2).*YHF(t1:t2)))/(sum(YHF(round(range)).*YHF(round(range))))*length(round(range))/nsta;
               %check if STA/LTA > threshold
               
        
%         
                if  ZHF(j)>=thHF; 
            if oth_HF ~= 1 %check if point before is over threshold
%                 fprintf(fid,'%6.6f %f\n',t1*dt, Z(j));
                Events_HF(cntrs,:)=YHF(t1:t2);
                
                Event_Time_HF(cntrs,:)=time(t1:t2);
                thresh_HF(cntrs)=ZHF(j);
                THF(j)=thHF;
                 S(:,:,cntrs)=data(t1-300:t1+1000,:);
                T=0:1e-6:(length(S)-1)*1e-6;
                T=T';
%                 [xn,yn]=locloc(S(:,:,cntrs),T);

               [xn,yn]=loci(S(:,:,cntrs),100,cntrs);
                LOC(:,cntrs)=[xn,yn];
                clear xn yn
                oth_HF=1;
               
                
                
                
                
                
                cntrs=cntrs+1;
            end
        else
            oth_HF=0 ;
        end
  
        
        
        
        
        
        
        
        
        
        

        
       
    end
    EvHF=Event_Time_HF(:,1);
    EV_HF=EvHF;

%   plot(t,ZHF)
%      hold on
%      plot(t,thHF,'r')
%      plot(t,THF,'r')
sig=S;
else
load detectedevents
S=sig;
dt=1e-6;
taperzone=0;
end

% % % % figure for a detected event
if evdetect
    for kk=1:23
        
intime(kk)=find(time>Event_Time_HF(kk,1),1);
    end
subplot(2,1,1)
set(gca,'fontsize',15)
osman=decimate(data(:,1),100);
tosman=decimate(time,100);
plot(tosman,osman,'linewidth',1.5)
hold on
for kk=1:23
plot([time(intime(kk)),time(intime(kk))],[-1 1],'r')
end
% plot([time(intime-300),time(intime-300)],[-1 1],'r')
% plot([time(intime+1000),time(intime+1000)],[-1 1],'r')
title('Detected Event Between Two Lines 0.0013 s')
xlabel('Time (s)')
ylabel('Amplitude (V)')

subplot(2,1,2)
plot(time(intime(3)-300)+T,S(:,1,3))
title('Zoomed Event')
xlabel('Time (s)')
ylabel('Amplitude (V)')
end

xrec=[-0.3 -0.3 -0.05 -0.05];
yrec=[0.05 -0.05 0.05 -0.05];
% ws=round(linspace(0,1,100)*700+50);
ws=round(linspace(0,1,50)*700+50);


if arrivcheck
thres=1;
for j=1:23
for si=1:length(xrec)
filtsig=taper(bp(sig(:,si,j),dt,3e4,5e4),taperzone/dt)  ; 
nsig=filtsig./max(abs(filtsig));
 AR=abs(nsig-mean(nsig))/std(nsig);
 AR(1:300)=0;
% AR=abs(filtsig-mean(filtsig))/std(filtsig);
if  sum(j==[4 7 9 10 11 13 15 16 17 19 20 21 23])
    AR=abs(filtsig-mean(filtsig))/std(filtsig);
AR(1:300)=0;
end
if sum(j==[ 8 22])
    AR=abs(filtsig-mean(filtsig))/std(filtsig);
    thres=0.35;
end
if sum(j==[ 18 ])
    AR=abs(filtsig-mean(filtsig))/std(filtsig);
    AR(1:100)=0;
    thres=1;
end
if sum(j==[ 6 ])
    AR=abs(filtsig-mean(filtsig))/std(filtsig);
    thres=0.4;
end
if sum(j==[ 2])
    AR=abs(filtsig-mean(filtsig))/std(filtsig);
    thres=0.55;
end
if sum(j==[1 ])
    AR=abs(filtsig-mean(filtsig))/std(filtsig);
    thres=0.05;
end
if sum(j==[3 ])
    AR=abs(filtsig-mean(filtsig))/std(filtsig);
    thres=0.25;
end
 if isempty(find(AR>thres,1))
     lim(si)=0;
     
 else
lim(si)=find(AR>thres,1);
% figure
% subplot(2,1,1)
% plot(sig(:,si,j))
% hold on
% plot(filtsig,'r')
% % hold all
% plot(lim(si)-20,filtsig(lim(si)-20),'g*')
% subplot(2,1,2)
% plot(AR)
% hold on
% plot(lim(si),AR(lim(si)),'g*')
% saveas(1,['arrtimethresh1_',num2str(j),'_',num2str(si)],'png')
% % pause
% close all
 end
end


Signal=sig(:,:,j);

presig=15;

        artime(j,1)=lim(1)+1;
if min(lim)<presig
    presig=0;
end

    
for i=1:length(ws)     
for wi=1:4
len=length(Signal(lim(wi)-presig:lim(wi)+1+ws(i),wi));
acc(1:len,wi)=Signal(lim(wi)-presig:lim(wi)+1+ws(i),wi);
end



    [res]=loci_v2(acc,Signal,lim);
%     pause
    clear acc
    xop(i,j)=res(1);
    yop(i,j)=res(2);
    xat(i,j)=res(3);
    yat(i,j)=res(4);
    xtr(i,j)=res(5);
    ytr(i,j)=res(6);
i


end
end
save('LocationRes.mat','xop','yop')



else
for j=1:23
for i=1:length(ws)
    [res]=loci_v2(sig(:,:,j),ws(i),0);
%     res=[xn;yn;xat;yat;xtr;ytr];
    xop(i,j)=res(1);
    yop(i,j)=res(2);
    xat(i,j)=res(3);
    yat(i,j)=res(4);
    xtr(i,j)=res(5);
    ytr(i,j)=res(6);
i
end

end


end





% % % % % % % % % refiltering and polarization check
polar=1;
if polar
for i=1:4
    for j=1:23
sfilt(:,i,j)=taper(bp(sig(:,i,j),1e-6,5e3,6e3),0);
    end
end
t=0:1e-6:1300*1e-6;
for j=1:size(sig,3)
for si=1:length(xrec)
    sfilt(:,si,j)=cumtrapz(cumtrapz(sfilt(:,si,j)));
% AR=abs(sfilt(:,si,j)-mean(sfilt(:,si,j)))/std(sfilt(:,si,j));
%  if isempty(find(AR>0.8,1))
%      lim(si)=0;
%      pol(si,j)=0;
%  else
% lim(si)=find(AR>0.8,1);
pol(si,j)=sfilt(800,si,j);
%  plot(sfilt(:,si,j))
%  hold all
%  plot(lim(si),sfilt(lim(si),si,j),'g*')
%  end

end
% saveas(1,['polcheck_',num2str(j)],'png')
%  close all
 npol(:,j)=pol(:,j)/max(abs(pol(:,j)));
end
end


% % % % % % % % %  end of polarization check

save('Uptohere.mat')
C = ColorBand( length(ws) );
Xrr=linspace(-0.4,0.35,903);
Yrr=linspace(0.15,-0.15,436);
t=0:1e-6:1e-6*1300;
for j=1:23
   
%     ccc=floor(EV_HF(j)*1e3);
%     cntr=floor(2+ccc/8);

    img=imread(['slide',num2str(j),'.png']);
    figure
  
  
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    subplot(2,1,1)
   imagesc(Xrr,Yrr,img(117:552,120:1022));

    %iptsetpref('ImshowAxesVisible','on');
    %imshow(img(117:552,120:1022), 'XData', Xrr, 'YData', Yrr)
    set(gca,'Ydir','normal')
    set(gca,'fontsize',22)
    colormap('gray')
%     colorbar
%     caxis([0 0.05])
    xlabel('Length (m)')
    ylabel('Length (m)')
    hold on
    scatter(xop(:,j),yop(:,j),60,C,'filled')

    
    for pp=1:4
        if  pol(pp,j)>=0
            POL(pp,:)=[1 0 1];
        else
            POL(pp,:)=[1 1 0];
        end
    end
    scatter(xrec(1),yrec(1),100,[1 0 1],'d','filled')
    scatter(xrec(2),yrec(2),100,[1 1 0],'d','filled')
    scatter(xrec,yrec,100,POL,'d','filled')

    
    legend('Source Location','Sensors pol (+)','Sensors pol (-)')
    axis([-0.4 0.35 -0.15 0.15])
%     subplot(3,1,2)
    
%     for k=1:length(ws)
%     inx=find(xop(k,j)<=Xrr,1);
%     iny=find(yop(k,j)<=Yrr,1);
%     immat=mean(im2double(img(117:552,120:1022)),3);
%     val(k)=(immat(iny,inx))/max(max(immat));
%     end
%     scatter(ws,val,20,C,'filled')
%     xlabel('Window Length (points dt=1E-6)')
%     ylabel('% Correspondance')
    subplot(2,1,2)
    plot(t,sig(:,1,j))
    set(gca,'fontsize',22)
    hold on
    plot(t(artime(j,1)+1),sig(artime(j,1)+1,1,j),'rv','markersize',10,'MarkerFaceColor','r')
    scatter(t(artime(j,1)+1+ws),sig(artime(j,1)+1+ws,1,j),40,C,'filled')
%     plot(t(artime(j,2)+1),sig(artime(j,2)+1,1,j),'g*')
    xlabel('Time (s)')
    ylabel('Amplitude (V)')
    legend('Signal','Arrival Time','Window End')
%     pause
% saveas(1,['loccheckhd_',num2str(j),'.eps'],'epsc')
saveas(1,['SuppMatEV_',num2str(j),'.eps'],'epsc')
saveas(1,['SuppMatEV_',num2str(j)],'png')
saveas(1,['SuppMatEV_',num2str(j)],'fig')
   close all
end
%    
%        plot(t,sig(:,1,j))
%     set(gca,'fontsize',18)
%     hold on
%     plot(t(artime(j,1)+1),sig(artime(j,1)+1,1,j),'r*')
%     scatter(t(artime(j,1)+1+ws),sig(artime(j,1)+1+ws,1,j),20,C,'filled')
% %     plot(t(artime(j,2)+1),sig(artime(j,2)+1,1,j),'g*')
%     xlabel('Time (s)')
%     ylabel('Amplitude (V)')
%     legend('Signal','Arrival Time','Window End')c
%  saveas(1,['3sub8Nov',num2str(j),'.eps'],'epsc')
% saveas(1,['3sub8Nov_',num2str(j)],'png')
%    close all  
%    figure
%    set(gca,'fontsize',18)
% xax = {'NW','SW','NE','SE'};
% bar(pol(:,j)/max(abs(pol(:,j))));
% set(gca, 'XTick', 1:4, 'XTickLabel', xax);
% title('Polarization check ("+" shows norm. compaction)')
%     saveas(1,['2sub8Nov',num2str(j),'.eps'],'epsc')
% saveas(1,['2sub8Nov_',num2str(j)],'png')
% close all
% end

% for j=1:23
% %     ccc=floor(EV_HF(j)*1e3);
% %     cntr=floor(2+ccc/8);
%     img=imread(['slide',num2str(j),'.png']);
%     figure
%   
%   
%     set(gcf,'units','normalized','outerposition',[0 0 1 1])
% %     subplot(3,1,1)
%     imagesc(Xrr,Yrr,img(117:552,120:1022));
%     set(gca,'Ydir','normal')
%     set(gca,'fontsize',18)
%     colormap('gray')
%     hold on
%     scatter(xop(:,j),yop(:,j),60,C,'filled')
% 
%     scatter(xat(:,j),yat(:,j),60,C,'s','filled')
%     scatter(xtr(:,j),ytr(:,j),60,C,'h','filled')
%     
%     for pp=1:4
%         if  pol(pp,j)>=0
%             POL(pp,:)=[1 0 1];
%         else
%             POL(pp,:)=[1 1 0];
%         end
%     end
%     scatter(xrec,yrec,60,POL,'d','filled')
% %         for ss=1:length(xop(:,j))
% %     for i=1:4
%     
% %     [xp,yp]=arrowLand(xop(ss,j),yop(ss,j),npol(i,j)/10,xrec(i),yrec(i));
% %     if npol(i,j)<0 
% %         xp=(-xop(ss,j)+xrec(i))*npol(i,j)/5;
% %         yp=(-yop(ss,j)+yrec(i))*npol(i,j)/5;
% %         quiver(xop(ss,j),yop(ss,j),xp,yp,'g')
% % %     else
% % %         xp=(-xop(ss,j)+xrec(i))*npol(i,j)/5;
% % %         yp=(-yop(ss,j)+yrec(i))*npol(i,j)/5;
% % %         xp1=xp+xop(ss,j);
% % %         yp1=yop(ss,j)+yp;
% % %         dx=-xp1+xop(ss,j);
% % %         dy=-yp1+yop(ss,j);
% % %         quiver(xp1,yp1,dx,dy,'--k')
% %     end
% %     if npol(i,j)>0
% %         xp=xop(ss,j)-xp;
% %         yp=yop(ss,j)-yp;
% %     
% %     else
% %         xp=xp-xop(ss,j);
% %         yp=yp-yop(ss,j);
% %     quiver(xp,yp,xop(ss,j),yop(ss,j),'--k')
% %     end
% %     end
% %     end
%     
%     legend('ESEH','ATL','TRL','Sensors')
%     axis([-0.4 0.35 -0.15 0.15])
% %     subplot(3,1,2)
% %     
% %     for k=1:length(ws)
% %     inx=find(xop(k,j)<=Xrr,1);
% %     iny=find(yop(k,j)<=Yrr,1);
% %     immat=mean(im2double(img(117:552,120:1022)),3);
% %     val(k)=(immat(iny,inx))/max(max(immat));
% %     end
% %     scatter(ws,val,20,C,'filled')
% %     xlabel('Window Length (points dt=1E-6)')
% %     ylabel('% Correspondance')
% % 
% 
% 
% %     pause
% % saveas(1,['loccheckhd_',num2str(j),'.eps'],'epsc')
% saveas(1,['4sub8Nov',num2str(j),'.eps'],'epsc')
% saveas(1,['4sub8Nov_',num2str(j)],'png')
% saveas(1,['4sub8Nov_',num2str(j)],'fig')
%    close all
% 
% end
% 
% imloc =' C:\Users\usr1\Documents\Dropbox Save 22012015\Experiments\December 13\Experiments 16122013\Movie\2 Bar Loose Packing';
% svloc = 'C:\Users\usr1\Dropbox\Experiments\June 15\Injection Exp';
% % svloc='C:\Users\usr1\Dropbox\Experiments\July 14\Traditional Beamforming\Dispersion_Reflection\Experimental_Signal\Injection Experiment\Project_Seer\7th Try Dispersion';
% imgname= '2exp_C001H001S0001000%03d.bmp';
% col = ColorBand(25);
% % col = linspace(1,10,1006);
% 
% xrec=[-0.3 -0.3 -0.05 -0.05];
% yrec=[-0.05 0.05 -0.05 0.05];
% Xr=(-0.4:0.005:0.4);
% Yr=(-0.2:0.005:0.2);
% window_length=1e-3;
% 
% con=1;
% for i=1:length(EV_HF)
%     ccc=floor(EV_HF(i)*1e3);
%     tic
%     cntr=floor(2+ccc/8);
%     img=imread(fullfile(imloc,sprintf(imgname, cntr)));
% 
% 
%     
%     figure(1)
%     subplot (2,1,1)
%     set(gcf,'units','normalized','outerposition',[0 0 1 1])
%     imagesc(Xr,Yr,img(107:645,1:1024));
%     colormap('gray')
%     hold on
% %     ex = exist(fullfile(svloc,[num2str(ccc),'_Exp.mat']));
% %     if ex;
% 
% %     load(fullfile(svloc,[num2str(i),'_Exp.mat']))
% %     locx(con)=xn;
% %     locy(con)=yn;
% %     con=con+1;
% %     end
% 
%     scatter(LOC(1,i),-LOC(2,i),'r')
% 
% 
% 
% 
%     plot(xrec,yrec,'go')
%     legend('Source','Sensors')
%     axis([-0.4 0.4 -0.2 0.2])
%    
%    title(['Acoustic Events During Experiment  T=',num2str(ccc*window_length),'seconds'])
%    subplot (2,1,2)
%    plot(T,S(:,1,con))
%    title('Signal')
%    savename=[num2str(con),'_deneme2'];
%    saveas(1,fullfile(svloc,savename),'png')
%    close all
%    con=con+1;
% %     mov(:,:,i)=getframe(i);
% %      close all
%      toc
%     end
% 
% 
% 
% 
% 

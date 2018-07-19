% clear all
% close all
% clc
% 
function [res]=loci_v2(acc,Tr_Acc,lim)
%  [err,res]=loci(xsrc,ysrc,acc,ws)
% % % % % % % % % % GLASS PLATE ONLY 
xrec=[-0.3 -0.3 -0.05 -0.05];
yrec=[0.05 -0.05 0.05 -0.05];
% function [err,res]=deneme_sim(xsrc,ysrc,reflections,broadband,dispersive)
% tic
% load localization06072015low_27
fig=0;
% 
% X=(-0.38:0.03:0.38);
% Y=(-0.19:0.05:0.19);
% addpath C:\Users\usr1\Dropbox\Super_Mega_Script_Folder
% addpath D:\Dropbox\Dropbox\Super_Mega_Script_Folder
addpath C:\Users\usr1\Dropbox\Super_Mega_Script_Folder
% addpath C:\Users\Dhareion\Dropbox\Super_Mega_Script_Folder
% xrec=[0 0.35 0.35 0];
% yrec=[-0.1 -0.1 0.1 0.1];
% load localization06072015low_better_27
% xsrc=0.10;
% ysrc=-0.09;
% % % Signal Generation Parameters
reflections=1;
broadband=0;
% winsize=1000;
Es_grid=10;
Xs=(-0.35:0.005:0.1);
Ys=(-0.15:0.005:0.15);
Fs=1e6;
% % % % % % % % % % GLASS PLATE ONLY 
t0=0:1/Fs:(length(acc(:,1))-1)*1/Fs;
rhob = 7800; %kg/m^3 densite 2500 verre / 7800 acier / 1140 polyamide / 2600 granite
Eb = 203e9; %Pa module d'Young 70e9 Verre / 203e9 acier / 4e9 polyamide / 60e9 granite
nub = 0.3; % coef Poisson 0.2 verre / 0.3 acier / 0.4 polyamide / 0.27 granite

rhop = 2500; %kg/m^3 2500 verre / 2800 marbre / 1180 gplexi
Ep = 70e9; %Pa 70e9 verre / 26e9 marbre / 2.4e9 gplexi
nup = 0.2; % 0.2 verre / 0.3 marbre / 0.375 gplexi

h = 0.01; %m % plate thickness
rayon = 0.001/2*2.5;% diametre de bille (mm)
m = 0.001*0.0638;% masse d'une bille d'acier (kg)
% m = 6.38e-6;
B = h^3*Ep/(12*(1-nup^2)); % bending stiffness m^3 Pa
E = 1/((1 - nub^2)/Eb + (1 - nup^2)/Ep); % (Pa) module equivalent
G=Ep/(2*(1+nup));  %shear modulus


% % % % % % % % % % GLASS PLATE ONLY 
% Signal=acc;
Tr_Acc=flipud(Tr_Acc);
% % % % % % % % % % % % % % % % % % %  Arrival Time Calc
% thres=1;
% for si=1:length(xrec)
% AR=abs(Signal(:,si)-mean(Signal(:,si)))/std(Signal(:,si));
%  if isempty(find(AR>thres,1))
%      lim(si)=0;
%      t_data(si)=0;
%  else
% lim(si)=find(AR>thres,1);
% t_data(si)=t0(lim(si));
% 
% 
%  end
% end
del=lim-min(lim);% % % % % % % % % % GLASS PLATE ONLY 



Loc=zeros(length(Ys),length(Xs));







% % % % % % % % % % % % % % % % % % %   Proposed Method 
   
presig=0;
% 
%     if min(lim)<presig
%         presig=min(lim)-10;
%     end
%     
% for i=1:4
% 
% len=length(Signal(lim(i)-presig:lim(i)+1+ws,i));
% acc(1:len,i)=Signal(lim(i)-presig:lim(i)+1+ws,i);
% 
% % plot(Signal(lim(i)-50:lim(i)+1+1000,i))
% % pause
% end
% Tr_Acc=flipud(Signal);
dt=1/Fs;
t0=0:dt:(length(acc(:,1))-1)*dt;



% % % % % % % % % % GLASS PLATE ONLY 
% imagesc(acc)
% pause
% hann1=hanning(length(t0));
for i=1:length(xrec)

[Fw,f0]=fftrl(acc(:,i),t0,0,2^nextpow2(length(t0)));
wave1(:,i)=Fw;
end
% omega = 2*pi*f0;
load dispexp
w1=omega;
k1=k;
omega = 2*pi*f0;
k=interp1(w1,k1,omega);
c=(diff(omega)./diff(k));
% plot(c)
% pause
% k = ((rhop*h/B)^(1/4).*sqrt(omega));% wave number
Vray=(0.87+1.12*nup)/(1+nup).*sqrt(G/rhop);
% Vph=omega./k;
% c=zeros(length(Vph));
% c=(diff(omega)./diff(k));
vind=find(c>Vray,1); %cut off frequency for non-dispersive phase
% % % % % % % % % % GLASS PLATE ONLY 
% c(vind:end)=Vray;
% k(vind:end)=k(vind)/omega(vind)*omega(vind:end);

meanfr=trapz(abs(Fw).*f0)./trapz(abs(Fw));
mind=find(f0>=meanfr,1);
fcutoff=f0(vind);
meanvel=c(vind);

% save('grvelocity.mat','c','f0')
% plot(gam)
% pause



% % % % % % % % % % GLASS PLATE ONLY 
  
        
 

% flim=find(f0>4e3,1);



for xi=1:length(Xs)
    for yi=1:length(Ys)
       
    ximg=Xs(xi);
    yimg=Ys(yi);

% % % % % % % % % % GLASS PLATE ONLY 
    %%%%%%%%%% Arrival Time Localization
    for pp=1:length(xrec)  
    R(pp)=sqrt((xrec(pp)-Xs(xi))^2+(yrec(pp)-Ys(yi))^2);
    tch(pp)=R(pp)/meanvel;
    end
    RMS(yi,xi)=sqrt(sum((((lim-min(lim))*1e-6)-(tch-min(tch))).^2));

    %%%%%%%%%%%%%%%%% Time Reversal
       for i=1:length(xrec)
       
            rimg=sqrt((xrec(i)-ximg)^2+(yrec(i)-yimg)^2);
            delay=rimg/meanvel;
            inf=round(delay*Fs)+1;
            sup=numel(Tr_Acc(:,i))+round(delay*Fs);
            highamp(inf:sup,i)=Tr_Acc(:,i);
       end
      
      Beam(yi,xi)=max(sum(highamp,2).^2);
%       plot(highamp)
%       pause
      clear highamp
      %%%%%%%%%%%%%%%%%%%%%
      
      
   % % % % % % % % % % GLASS PLATE ONLY    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Energy Proposed

        for ii=1:length(xrec)




                rimg=sqrt((xrec(ii)-ximg)^2+(yrec(ii)-yimg)^2);
                
                
                 
                TRAP=abs(wave1(:,ii)).^2./omega(:).^2;
                TRAP(1)=0;
                Eimg(ii)=rhop*h*trapz(TRAP(3:end).*c(1:end-1));
                Ech(ii)=rimg*Eimg(ii);
 
        end
        NormE=Ech./max(Ech);
            Loc_new(yi,xi)=std(NormE);
%            Loc_new(yi,xi)=sqrt(1/6*((Eimg(1)-Eimg(2))^2+(Eimg(1)-Eimg(3))^2+(Eimg(1)-Eimg(4))^2+(Eimg(2)-Eimg(3))^2+...
%            (Eimg(2)-Eimg(4))^2+(Eimg(3)-Eimg(4))^2));
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %%%%%%%%%%%%%%%%% Energy Conventional

   % % % % % % % % % % GLASS PLATE ONLY 
         
end
    
end
% figure
% imagesc(Loc_new)
% pause
%%% Energy Prop
[j1,j2]=find(Loc_new==min(min(Loc_new)));
xn=Xs(j2(1));
yn=Ys(j1(1));
% im1=Loc_new./max(max(abs(Loc_new)));

% %%%% Time Reversal
Nbeam=Beam./max(max(Beam));
[tr1,tr2]=find(Nbeam==max(max(Nbeam)));
% imagesc(Nbeam)
% pause
xtr=Xs(tr2(1));
ytr=Ys(tr1(1));
% im3=1./Nbeam;
% im3=im3./max(max(im3));
% %%%% Arrival Time
[at1,at2]=find(RMS==min(min(RMS)));
xat=Xs(at2(1));
yat=Ys(at1(1));
% im4=RMS./max(max(abs(RMS)));

% 
% econ=sqrt((xc-xsrc)^2+(yc-ysrc)^2);
% eat=sqrt((xat-xsrc)^2+(yat-ysrc)^2);
% etr=sqrt((xtr-xsrc)^2+(ytr-ysrc)^2);
% en=sqrt((xn-xsrc)^2+(yn-ysrc)^2);
% for i=1:length(xn)
% en(i)=sqrt((xn(i)-xsrc)^2+(yn(i)-ysrc)^2);
% end
% eni=min(en);
% eno=find(en==min(en));
% THIS CODE HAS PARAMETERS ONLY FOR PLEXIGLASS PLATE !!!!!!!!

% err(j)=en;
% round(j/30*100)


% err=[en;eat;econ;etr];
if fig
    figure
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
 subplot(2,2,1)
    imagesc(Xs,Ys,im1)
    colorbar
    hold on
    
    plot(xrec,yrec,'r.','MarkerSize',20)
%     plot(xsrc,ysrc,'m.','MarkerSize',20)
    plot(xn,yn,'g.','MarkerSize',20)
%     caxis([0 1])
      title(['PEBL Error='])
    subplot(2,2,2)
    imagesc(Xs,Ys,im4)
    colorbar

    hold on
    plot(xrec,yrec,'r.','MarkerSize',20)
%     plot(xsrc,ysrc,'m.','MarkerSize',20)
    plot(xat,yat,'g.','MarkerSize',20)
%     caxis([0 1])
    title(['ATL Error='])

    subplot(2,2,4) 
    imagesc(Xs,Ys,im3)
    colorbar
  
     hold on
    
    plot(xrec,yrec,'r.','MarkerSize',20)
%     plot(xsrc,ysrc,'m.','MarkerSize',20)
    plot(xtr,ytr,'g.','MarkerSize',20)
%     caxis([0 1])
    title(['TRL Error='])
%     sname=['Detection_',num2str(fig),'.png'];
% export_fig(gcf,sname)
%     saveas(1,['detect_',num2str(fig)],'png')
    pause
    close all
end
% figure
% plot(t0,acc)
res=[xn;yn;xat;yat;xtr;ytr];
end
% toc
% toc
% for i=1:4
%     RR(i)=sqrt((xs(27)-xrec(i))^2+(yrec(i)-ys(27))^2);
% end
% figure
% subplot(2,2,1)
% imagesc(squeeze(Acc(:,:,1)))
% colorbar
% axis([50 90 50 400])
% caxis([-0.2 0.2])
% subplot(2,2,2)
% imagesc(squeeze(Acc(:,:,2)))
% colorbar
% axis([50 90 50 400])
% caxis([-0.2 0.2])
% subplot(2,2,3)
% imagesc(squeeze(Acc(:,:,3)))
% colorbar
% axis([50 90 50 400])
% caxis([-0.2 0.2])
% subplot(2,2,4)
% imagesc(squeeze(Acc(:,:,4)))
% colorbar
% caxis([-0.2 0.2])
% axis([50 90 50 400])
 %%%% This program reads h3dtdinv/h3doctree observations and predicted data files and
%%%% generates a profile plot. Only valid for observation data located on 
%%%% a single line i.e for 2D inversion using the 3D UBC EM codes  

%%%% Update June 2014. Ability to detect restarted inversions in an INV2 directory 

%% open files




clear all;


workpath='C:\LC\Private\dominiquef\Projects\4234_Anglo_Mosku_3DTM\Block14_Inv1\';
workfile=fopen([workpath,'Work_list.dat']);

datapath='C:\LC\Private\dominiquef\Projects\4234_Anglo_Mosku_3DTM\Block14_Figures'

% Line oriented EW = 0 , NS = 1
flag = 0;

%loop over directories

while 1
close all;    
line=fgetl(workfile);

if line==-1
     break; break;
else
line=strtrim(line);    
obspathname=[workpath,line];

eval(['cd ',obspathname]);

dir_contents=ls;
S=size(dir_contents);


for ii =1:1:S(1)
    if strmatch(strtrim(dir_contents(ii,:)),'INV2')
        cd INV2
        dir_contents=ls;
        obspathname=pwd;
        %keyboard;
        break
    end
end



%find .obs file

sw_findobs=0;
S=size(dir_contents);
for ii=3:1:S(1)
    s=strtrim(dir_contents(ii,:));
    if strcmp(s(end-3:end),'.obs')
        sw_findobs=1;
        obs=s;
    end
end

if sw_findobs==0
    fprintf('\nFailed to locate file with .obs extension in inversion directory: \n %s',obspathname); 
    keyboard
end

    
        

%find most recent pred file

pred_num=0;    
 for ii=3:1:S(1)  
     s=strtrim(dir_contents(ii,:));
     if length(s) < 5 
     else
     if strcmp(s(1:5),'dpred')==1
         num=str2num(s(7:8));
         if pred_num < num
             pred_num=num;
             pre=s;
         end
     end
     end
 end
 

 if pred_num==0
    %keyboard
     fprintf('\n\n\nFailed to locate any predicted data files in inversion directory: \n %s',obspathname); 
    continue
    
 end
 
 
% find the out file

for ii=3:1:S(1)
    s=strtrim(dir_contents(ii,:));
    if strcmp(s(end-3:end),'.out')==1
        out=s;
    end
end

 
 
name=[obspathname,'\',out]    
    
A=load(name);
SA=size(A)

if SA(1) > 0
scrsz = get(0,'ScreenSize');
figure('Position',[1 1 scrsz(3) scrsz(4)])


subplot(2,2,1)
plot(1:1:SA(1), A(:,3),'.b')
grid on
title('Data Misfit vs iteration - Linear')
xlabel('iteration');
ylabel('data misfit');

subplot(2,2,2)
semilogy(1:1:SA(1), A(:,3),'.b')
grid on
title('Data Misfit vs iteration - Log')
xlabel('iteration');
ylabel('data misfit');

subplot(2,2,3)
semilogx(A(:,1),A(:,3),'.b')
grid on
title('Data Misfit vs tradeoff parameter')
xlabel('Beta');
ylabel('Data Misfit');

subplot(2,2,4)
loglog(A(:,1),A(:,3),'.b')
grid on
title('Data Misfit vs tradeoff parameter')
xlabel('Beta');
ylabel('Data Misfit');

h=gcf;
s=[obspathname,'\',obs,'_data_misfit.jpg'];
print( h, '-djpeg90', s);

h=gcf;
s=[datapath,'\',obs,'_data_misfit.jpg'];
print( h, '-djpeg90', s);

end

%% Open the files
fido=fopen([obspathname,'\',obs],'r');
fidp=fopen([obspathname,'\',pre],'r');



%% read data

clear obsdata
clear predata

stationcount=0;
while 1
    linecount=0;
    tline = fgetl(fido);
    if isempty(tline)
    elseif ~ischar(tline)
        break;
    elseif strcmp(tline(1:6),'N_TIME')
        ntime=str2num(tline(7:end));
        stationcount=stationcount+1;
        while linecount < ntime
            tline = fgetl(fido);



            if isempty(tline)
            else
                linecount=linecount+1;
                obsdata(stationcount,linecount,:)=str2num(tline);
                
                %read a line from the pred file

                tlinep = fgetl(fidp);
                if isempty(tlinep)
                    while isempty(tlinep)
                        tlinep = fgetl(fidp);
                    end
                end
                if ischar(tlinep)
                    predata(stationcount,linecount,:)=str2num(tlinep);
                end

            end
        end
    end
end


%% convert units
obsdata(:,:,4:end)=obsdata(:,:,4:end)*1e12; % convert to pT per A
predata(:,:,4:end)=predata(:,:,4:end)*1e12; % convert to pT per A


fclose(fido); 
fclose(fidp);

%% plot Z
dropchannel=13;
%dropchannel=dropchannel-2;
S=size(predata)
scrsz = get(0,'ScreenSize');
figure('Position',[1 1 scrsz(3) scrsz(4)])

ntime_perpannel=floor(S(2)/3);

extratime=abs(ntime_perpannel*3 - S(2));

subplot(3,1,1);
for ii=1:1:ntime_perpannel+extratime
    plot(obsdata(:,ii,2-flag),obsdata(:,ii,21),'k')
   % plot(obsdata(:,ii,2),obsdata(:,ii,21),'.k')
    hold on
    plot(predata(:,ii,2-flag),predata(:,ii,13),'r')
end
xlabel('northing [m]');
ylabel('pV m^{2} / A');
titlestr=['Z comp Channel ',num2str(dropchannel+2),' to ',num2str(2*(ntime_perpannel+extratime)+dropchannel)];
title(titlestr);
grid on
grid minor


subplot(3,1,2)
for ii=ntime_perpannel+extratime:1:ntime_perpannel*2+extratime
    plot(obsdata(:,ii,2-flag),obsdata(:,ii,21),'k')
    hold on
    %plot(obsdata(:,ii,2),obsdata(:,ii,21),'.k')
    plot(predata(:,ii,2-flag),predata(:,ii,13),'r')
end
xlabel('northing [m]');
ylabel('pV m^{2} / A');
titlestr=['Z comp Channel ',num2str(2*(ntime_perpannel+extratime)+dropchannel), ' to ',num2str(ntime_perpannel*4+dropchannel)];
title(titlestr);
grid on
grid minor

subplot(3,1,3)
for ii=ntime_perpannel*2+extratime:1:S(2)
    plot(obsdata(:,ii,2-flag),obsdata(:,ii,21),'k')
    hold on
    %plot(obsdata(:,ii,2),obsdata(:,ii,21),'.k')
    plot(predata(:,ii,2-flag),predata(:,ii,13),'r')
end
xlabel('northing [m]');
ylabel('pV m^{2} / A');
titlestr=['Z comp Channel ',num2str(ntime_perpannel*4+dropchannel), ' to ',num2str(S(2)*2+dropchannel)];
title(titlestr);
grid on
grid minor

h=gcf;
s=[obspathname,'\',obs,'_Z_comp.jpg'];
print( h, '-djpeg90', s);

h=gcf;
s=[datapath,'\',obs,'_Z_comp.jpg'];
print( h, '-djpeg90', s);

%%
% %% plot X
% figure(2)
% subplot(3,1,1);
% for ii=1:1:ntime_perpannel+extratime
%     plot(obsdata(:,ii,2),obsdata(:,ii,17),'k')
%     hold on
%     plot(predata(:,ii,2),predata(:,ii,11),'r')
% end
% xlabel('northing [m]');
% ylabel('pV m^{2} / A');
% titlestr=['X comp Channel ',num2str(dropchannel+2),' to ',num2str(2*(ntime_perpannel+extratime)+dropchannel)];
% title(titlestr);
% grid on
% grid minor
% 
% subplot(3,1,2)
% for ii=ntime_perpannel+extratime:1:ntime_perpannel*2+extratime
%     plot(obsdata(:,ii,2),obsdata(:,ii,17),'k')
%     hold on
%     plot(predata(:,ii,2),predata(:,ii,11),'r')
% end
% xlabel('northing [m]');
% ylabel('pV m^{2} / A');
% titlestr=['X comp Channel ',num2str(ntime_perpannel+extratime), ' to ',num2str(ntime_perpannel*2)];
% title(titlestr);
% grid on
% grid minor
% 
% subplot(3,1,3)
% for ii=ntime_perpannel*2+extratime:1:S(2)
%     plot(obsdata(:,ii,2),obsdata(:,ii,17),'k')
%     hold on
%     plot(predata(:,ii,2),predata(:,ii,11),'r')
% end
% xlabel('northing [m]');
% ylabel('pV m^{2} / A');
% titlestr=['X comp Channel ',num2str(ntime_perpannel*4+dropchannel), ' to ',num2str(S(2)*2 +dropchannel)];
% title(titlestr);
% grid on
% grid minor
% 
% 
% h=gcf;
% s=[obspathname,obs,'_X_comp.jpg']
% print( h, '-djpeg90', s);
% 


%% plot Y
scrsz = get(0,'ScreenSize');
figure('Position',[1 1 scrsz(3) scrsz(4)])
subplot(3,1,1);
for ii=1:1:ntime_perpannel+extratime
    plot(obsdata(:,ii,2-flag),obsdata(:,ii,19-2*flag),'k')
    hold on
    plot(predata(:,ii,2-flag),predata(:,ii,12-flag),'r')
end
xlabel('northing [m]');
ylabel('pV m^{2} / A');
titlestr=['Inline comp Channel ',num2str(dropchannel+2),' to ',num2str(2*(ntime_perpannel+extratime)+dropchannel)];
title(titlestr);
grid on
grid minor

subplot(3,1,2)
for ii=ntime_perpannel+extratime:1:ntime_perpannel*2+extratime
    plot(obsdata(:,ii,2-flag),obsdata(:,ii,19-2*flag),'k')
    hold on
    plot(predata(:,ii,2-flag),predata(:,ii,12-flag),'r')
end
xlabel('northing [m]');
ylabel('pV m^{2} / A');
titlestr=['Inline comp Channel ',num2str(2*(ntime_perpannel+extratime)+dropchannel), ' to ',num2str(ntime_perpannel*4+dropchannel)];
title(titlestr);
grid on
grid minor

subplot(3,1,3)
for ii=ntime_perpannel*2+extratime:1:S(2)
    plot(obsdata(:,ii,2-flag),obsdata(:,ii,19-2*flag),'k')
    hold on
    plot(predata(:,ii,2-flag),predata(:,ii,12-flag),'r')
end
xlabel('northing [m]');
ylabel('pV m^{2} / A');
titlestr=['Inline comp Channel ',num2str(ntime_perpannel*4+dropchannel), ' to ',num2str(S(2)*2+dropchannel)];
title(titlestr);
grid on
grid minor


h=gcf;
s=[obspathname,'\',obs,'_Inline_comp.jpg'];
print( h, '-djpeg90', s);

h=gcf;
s=[datapath,'\',obs,'_Inline_comp.jpg'];
print( h, '-djpeg90', s);


%% plot Z
dropchannel=13;
%dropchannel=dropchannel-2;
S=size(predata)
scrsz = get(0,'ScreenSize');
figure('Position',[1 1 scrsz(3) scrsz(4)])

ntime_perpannel=floor(S(2)/3);

extratime=abs(ntime_perpannel*3 - S(2));

subplot(3,1,1);
for ii=1:1:ntime_perpannel+extratime
    semilogy(obsdata(:,ii,2-flag),obsdata(:,ii,21),'k')
   % plot(obsdata(:,ii,2),obsdata(:,ii,21),'.k')
    hold on
    semilogy(predata(:,ii,2-flag),predata(:,ii,13),'r')
end
xlabel('northing [m]');
ylabel('pV m^{2} / A');
titlestr=['Z comp Channel ',num2str(dropchannel+2),' to ',num2str(2*(ntime_perpannel+extratime)+dropchannel)];
title(titlestr);
grid on
grid minor


subplot(3,1,2)
for ii=ntime_perpannel+extratime:1:ntime_perpannel*2+extratime
    semilogy(obsdata(:,ii,2-flag),obsdata(:,ii,21),'k')
    hold on
    %plot(obsdata(:,ii,2),obsdata(:,ii,21),'.k')
    semilogy(predata(:,ii,2-flag),predata(:,ii,13),'r')
end
xlabel('northing [m]');
ylabel('pV m^{2} / A');
titlestr=['Z comp Channel ',num2str(2*(ntime_perpannel+extratime)+dropchannel), ' to ',num2str(ntime_perpannel*4+dropchannel)];
title(titlestr);
grid on
grid minor

subplot(3,1,3)
for ii=ntime_perpannel*2+extratime:1:S(2)
    semilogy(obsdata(:,ii,2-flag),obsdata(:,ii,21),'k')
    hold on
    %plot(obsdata(:,ii,2),obsdata(:,ii,21),'.k')
    semilogy(predata(:,ii,2-flag),predata(:,ii,13),'r')
end
xlabel('northing [m]');
ylabel('pV m^{2} / A');
titlestr=['Z comp Channel ',num2str(ntime_perpannel*4+dropchannel), ' to ',num2str(S(2)*2+dropchannel)];
title(titlestr);
grid on
grid minor

h=gcf;
s=[obspathname,'\',obs,'_Z_comp_log_scale.jpg'];
print( h, '-djpeg90', s);

h=gcf;
s=[datapath,'\',obs,'_Z_comp_log_scale.jpg'];
print( h, '-djpeg90', s);



%% plot Y
scrsz = get(0,'ScreenSize');
figure('Position',[1 1 scrsz(3) scrsz(4)])
subplot(3,1,1);
for ii=1:1:ntime_perpannel+extratime
    plot(obsdata(:,ii,2-flag),asinh(obsdata(:,ii,19-2*flag)),'k')
    hold on
    plot(predata(:,ii,2-flag),asinh(predata(:,ii,12-flag)),'r')
end
xlabel('northing [m]');
ylabel('arcsinh pV m^{2} / A');
titlestr=['Inline comp Channel ',num2str(dropchannel+2),' to ',num2str(2*(ntime_perpannel+extratime)+dropchannel)];
title(titlestr);
grid on
grid minor

subplot(3,1,2)
for ii=ntime_perpannel+extratime:1:ntime_perpannel*2+extratime
    plot(obsdata(:,ii,2-flag),asinh(obsdata(:,ii,19-2*flag)),'k')
    hold on
    plot(predata(:,ii,2-flag),asinh(predata(:,ii,12-flag)),'r')
end
xlabel('northing [m]');
ylabel('arcsinh pV m^{2} / A');
titlestr=['Inline comp Channel ',num2str(2*(ntime_perpannel+extratime)+dropchannel), ' to ',num2str(ntime_perpannel*4+dropchannel)];
title(titlestr);
grid on
grid minor

subplot(3,1,3)
for ii=ntime_perpannel*2+extratime:1:S(2)
    plot(obsdata(:,ii,2-flag),asinh(obsdata(:,ii,19-2*flag)),'k')
    hold on
    plot(predata(:,ii,2-flag),asinh(predata(:,ii,12-flag)),'r')
end
xlabel('northing [m]');
ylabel('arcsinh pV m^{2} / A');
titlestr=['Inline comp Channel ',num2str(ntime_perpannel*4+dropchannel), ' to ',num2str(S(2)*2+dropchannel)];
title(titlestr);
grid on
grid minor 


h=gcf;
s=[obspathname,'\',obs,'_Inline_comp_arcsinh_scale.jpg'];
print( h, '-djpeg90', s);

h=gcf;
s=[datapath,'\',obs,'_Inline_comp_arcsinh_scale.jpg'];
print( h, '-djpeg90', s);

end
end
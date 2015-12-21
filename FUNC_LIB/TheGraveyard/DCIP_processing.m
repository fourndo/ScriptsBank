clear all;
close all;

cd C:\PVK_Projects\3811_Copperfox\Processing\data_raw

survey=load('L400N-Run2.txt');

t1=survey(:,1);
r1=survey(:,2);
r2=survey(:,3);

for i=1:length(r2)
vp(i)=survey(i,5);
end
for i=1:length(r2)
ip(i)=survey(i,6)/1000;
end
s=size(vp);
si=size(ip);
vp_res=reshape(vp,s(2),1);
ip_res=reshape(ip,s(2),1);

cd C:\PVK_Projects\3811_Copperfox\Processing\data_edited

fob=fopen('dc_L400N-Run2_local_forward.dat','w');
fprintf(fob,'UBC-GIF DC Model\npole-dipole');
 for q=1:length(r2)
     if ((t1(q) < r1(q)) && r1(q) ~= r2(q))
    fprintf(fob,'\n\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f ',t1(q), t1(q), r1(q), r2(q), (vp_res(q)));
     end
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
fob=fopen('dc_L400N-Run2_local_all.dat','w');
fprintf(fob,'UBC-GIF DC Model\npole-dipole');
 for q=1:length(r2)
     if (((t1(q) > r2(q)) && (r1(q) ~= r2(q)) || (t1(q) < r1(q)) && (r1(q) ~= r2(q))))
    fprintf(fob,'\n\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f ',t1(q), t1(q), r1(q), r2(q), (vp_res(q)));
     end
 end
fclose('all');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fob=fopen('dc_L400N-Run2_local_back.dat','w');
fprintf(fob,'UBC-GIF DC Model\npole-dipole');
 for q=1:length(r2)
     if (((t1(q) > r2(q)) && r1(q) ~= r2(q)))
    fprintf(fob,'\n\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f ',t1(q), t1(q), r1(q), r2(q), (vp_res(q)));
     end
 end
fclose('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fob=fopen('ip_L400N-Run2_local_back.dat','w');
fprintf(fob,'UBC-GIF IP Model\npole-dipole');
 for q=1:length(r2)
     if ((t1(q) > r2(q) && r1(q) ~= r2(q)))
        fprintf(fob,'\n\t%8.11f\t%8.11f\t%8.11f\t%8.11f\t%8.11f ',t1(q), t1(q), r1(q), r2(q), (ip_res(q)));
     end
 end
fclose('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fob=fopen('ip_L400N-Run2_local_forward.dat','w');
fprintf(fob,'UBC-GIF IP Model\npole-dipole');
 for q=1:length(r2)
     if ((t1(q) < r1(q) && r1(q) ~= r2(q)))
        fprintf(fob,'\n\t%8.11f\t%8.11f\t%8.11f\t%8.11f\t%8.11f ',t1(q), t1(q), r1(q), r2(q), (ip_res(q)));
     end
 end
fclose('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fob=fopen('ip_L400N-Run2_local_all.dat','w');
fprintf(fob,'UBC-GIF IP Model\npole-dipole');
for q=1:length(r2)
     if (((t1(q) > r2(q) && r1(q) ~= r2(q))) || (t1(q) < r1(q) && r1(q) ~= r2(q)) )
    fprintf(fob,'\n\t%8.11f\t%8.11f\t%8.11f\t%8.11f\t%8.11f ',t1(q), t1(q), r1(q), r2(q), (ip_res(q)));
     end
 end
fclose('all');



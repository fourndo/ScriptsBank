% Program DCIP2UBC
% This program generates DCIP2D lines in UBC and Geosoft format from the 
% raw data files (Quantec format). The lines are seperated in the forward
% and backward survey (Titan lines). A prefix is added to the file name to
% differentiate between IP and DC lines. This prefix is used by the
% following scripts:
%
% DC_Batch_Inv.m : Batch file to invert the DC files for conductivity
%
% IP_Batch_Inv.m : Batch file to invert the IP files for chargeability
%
% DCIP_Batch_STD.m: Batch file to automatically assign a standard deviation
% to the data from the UBC log file after inversion.
%
% DCIP_Batch_DOI.m: Batch file to invert for 2 different reference models
% that can be used for the Depth of Investigation (DOI) calculations.
%
%
% INPUT:
% CSV files to be processed in a folder.
% Path of the INPUT folder containing all raw files in CSV format.
% Path for the output folder.

clear all;
close all;

home_dir = pwd;
%% User INPUT>>>>>
input_folder = 'C:\Projects\4001_Wallbridge_SkynnerLake_DCIP\Processing\Raw_files';
output_folder = 'C:\Projects\4001_Wallbridge_SkynnerLake_DCIP\Processing\Extracted_lines';

%% Driver

file_list=ls(input_folder);
cd(input_folder)

nb_files=size(file_list,1)-2;

% Cycle through all the file in the INPUT folder
for oo=1:nb_files
    
    AA= strtok(file_list(oo+2,:),' ');

    datafile=importdata(AA);

    survey=datafile.data;
    t1=survey(:,1);
    r1=survey(:,5);
    r2=survey(:,7);
    
    % Scale potential on current for UBC format
        for ii=1:length(r2)
            vp(ii)=survey(ii,11);
        end
    
    % Convert phaze from milli-rad -> rad   
        for i=1:length(r2)
            ip(i)=survey(i,13)/1000;
        end
        
    % Reshape vectors
%     s=size(vp);
%     si=size(ip);
%     vp_res=reshape(vp,s(2),1);
%     ip_res=reshape(ip,s(2),1);
    vp = vp(:);
    ip = ip(:);
    
    % Take the name of the file without the extension
    linename=AA(1:end-4);

    fob=fopen([output_folder '\dc_' linename '_forward.dat'],'w');
    fprintf(fob,'UBC-GIF DC Model\npole-dipole');
         for q=1:length(r2)
             if ((t1(q) < r1(q)) && r1(q) ~= r2(q))
            fprintf(fob,'\n\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f ',t1(q), t1(q), r1(q), r2(q), (vp(q)));
             end
         end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% fob=fopen([output_folder '\dc_' linename '_all.dat'],'w');
% fprintf(fob,'UBC-GIF DC Model\npole-dipole');
%  for q=1:length(r2)
%      if (((t1(q) > r2(q)) && (r1(q) ~= r2(q)) || (t1(q) < r1(q)) && (r1(q) ~= r2(q))))
%     fprintf(fob,'\n\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f ',t1(q), t1(q), r1(q), r2(q), (vp(q)));
%      end
%  end
% fclose('all');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fob=fopen([output_folder '\dc_' linename '_back.dat'],'w');
fprintf(fob,'UBC-GIF DC Model\npole-dipole');
 for q=1:length(r2)
     if (((t1(q) > r2(q)) && r1(q) ~= r2(q)))
    fprintf(fob,'\n\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f ',t1(q), t1(q), r1(q), r2(q), (vp(q)));
     end
 end
fclose('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fob=fopen([output_folder '\ip_' linename '_back.dat'],'w');
fprintf(fob,'UBC-GIF IP Model\npole-dipole');
 for q=1:length(r2)
     if ((t1(q) > r2(q) && r1(q) ~= r2(q)))
        fprintf(fob,'\n\t%8.11f\t%8.11f\t%8.11f\t%8.11f\t%8.11f ',t1(q), t1(q), r1(q), r2(q), (ip(q)));
     end
 end
fclose('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fob=fopen([output_folder '\ip_' linename '_forward.dat'],'w');
fprintf(fob,'UBC-GIF IP Model\npole-dipole');
 for q=1:length(r2)
     if ((t1(q) < r1(q) && r1(q) ~= r2(q)))
        fprintf(fob,'\n\t%8.11f\t%8.11f\t%8.11f\t%8.11f\t%8.11f ',t1(q), t1(q), r1(q), r2(q), (ip(q)));
     end
 end
fclose('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fob=fopen([output_folder '\Geosoft_' linename '_all.dat'],'w');
fprintf(fob,'Titan-Geosoft Format\n');
fprintf(fob,['LINE:' AA(1:(end-10)) ' DIPOLE:100 UNITS:M T=540,1370\n']);
fprintf(fob,'C1-x	C1-y	C2-x	C2-y	P1-x	P1-y	P2-x	P2-y	Current	CurrentErr	Vp	VpErr	Phz	PhzErr	AppRes');
for qq=1:size(datafile.data,1)
    fprintf(fob,'\n');
    for ii=1:size(datafile.data,2)-1
    fprintf(fob,'%8.11f,',datafile.data(qq,ii));
    end
    fprintf(fob,'%8.11f',datafile.data(qq,end));
 end
fclose('all');

end

cd(home_dir);
function [dhid,dobs,dpred,dstd,dres,dloc] = H3D_Pred_vs_obs(Atlas,work_dir,obsfile,predfile,index,argin)
% Read obs and predicted data and write to file 
% INPUT
% Data_DH: Structure containing HoleID.id and HoleID.depth for every
% observation
%           
% work_dir: directory where to find data.obs and data.pre
% obsfile: name of observation data file
% predfile: name of predicted data file
%
% OUTPUT
% HOLEid, Depth, data1, data2, ...
% Repeat for dBx, dBy, dBz


rfid1=fopen([work_dir '\' obsfile],'rt');
rfid2=fopen([work_dir '\' predfile],'rt');

if strcmp(argin,'write')==1
    fprintf('files will be written in Gocad format, might take a while\n');
    wobsx= fopen([work_dir '\' 'DHEM_dBx.obs'],'w');
    wobsy= fopen([work_dir '\' 'DHEM_dBy.obs'],'w');
    wobsz= fopen([work_dir '\' 'DHEM_dBz.obs'],'w');

    wpredx= fopen([work_dir '\' 'DHEM_dBx.pre'],'w');
    wpredy= fopen([work_dir '\' 'DHEM_dBy.pre'],'w');
    wpredz= fopen([work_dir '\' 'DHEM_dBz.pre'],'w');

    wresx= fopen([work_dir '\' 'DHEM_dBx.res'],'w');
    wresy= fopen([work_dir '\' 'DHEM_dBy.res'],'w');
    wresz= fopen([work_dir '\' 'DHEM_dBz.res'],'w');
end


line=fgets(rfid1);
new_id = 0; 

while isempty(strfind(line,'N_RECV'))==1
                line=fgets(rfid1);
end

nstn = str2num(line(7:end));
line=fgets(rfid1);
ntime = str2num(line(7:end));

dhid = zeros(nstn,1);

dpred{1} = zeros(nstn,ntime);
dpred{2} = zeros(nstn,ntime);
dpred{3} = zeros(nstn,ntime);

dobs{1} = zeros(nstn,ntime);
dobs{2} = zeros(nstn,ntime);
dobs{3} = zeros(nstn,ntime);

dstd{1} = zeros(nstn,ntime);
dstd{2} = zeros(nstn,ntime);
dstd{3} = zeros(nstn,ntime);

dres{1} = zeros(nstn,ntime);
dres{2} = zeros(nstn,ntime);
dres{3} = zeros(nstn,ntime);
 
dloc = zeros(nstn,3);

for ii = 1 : nstn
    
        dhid(ii) = Atlas.Htag(index(ii));
        

    for jj = 1 : ntime
        
        line1 = str2num(fgets(rfid1));
        line2 = str2num(fgets(rfid2));
        
        while length(line1)<2
            
            line1 = str2num(fgets(rfid1));
            
        end
        
        while length(line2)<2
            
            line2 = str2num(fgets(rfid2));
            
        end
        
        dloc(ii,1:3) = line1(1:3);
        
        dpred{1}(ii,jj) = line2(11);
        dpred{2}(ii,jj) = line2(12);
        dpred{3}(ii,jj) = line2(13);

        dobs{1}(ii,jj) = line1(17);
        dobs{2}(ii,jj) = line1(19);
        dobs{3}(ii,jj) = line1(21);

        dstd{1}(ii,jj) = line1(18);
        dstd{2}(ii,jj) = line1(20);
        dstd{3}(ii,jj) = line1(22);
        
        dres{1}(ii,jj) = abs( dobs{1}(ii,jj)-dpred{1}(ii,jj) ) / dstd{1}(ii,jj);
        dres{2}(ii,jj) = abs( dobs{2}(ii,jj)-dpred{2}(ii,jj) ) / dstd{2}(ii,jj);
        dres{3}(ii,jj) = abs( dobs{3}(ii,jj)-dpred{3}(ii,jj) ) / dstd{3}(ii,jj);
        
        if strcmp(argin,'write')==1
        %% Write data to file for re-survey
        % Add heather first
        if ii == 1 && jj == 1

            fprintf(wobsx,'HID DEPTH ');
            fprintf(wobsy,'HID DEPTH ');
            fprintf(wobsz,'HID DEPTH ');

            fprintf(wpredx,'HID DEPTH ');
            fprintf(wpredy,'HID DEPTH ');
            fprintf(wpredz,'HID DEPTH ');

            fprintf(wresx,'HID DEPTH ');
            fprintf(wresy,'HID DEPTH ');
            fprintf(wresz,'HID DEPTH ');

            for tt = 1:ntime

                fprintf(wobsx,'Time%i ',tt);
                fprintf(wobsy,'Time%i ',tt);
                fprintf(wobsz,'Time%i ',tt);

                fprintf(wpredx,'Time%i ',tt);
                fprintf(wpredy,'Time%i ',tt);
                fprintf(wpredz,'Time%i ',tt);

                fprintf(wresx,'Time%i ',tt);
                fprintf(wresy,'Time%i ',tt);
                fprintf(wresz,'Time%i ',tt);
                
            end

        end
        
        if jj ==1
            
            fprintf(wobsx,'\n%s %8.3f ',char(Atlas.id{index(ii)}),Atlas.DXYZ(index(ii)));
            fprintf(wobsy,'\n%s %8.3f ',char(Atlas.id{index(ii)}),Atlas.DXYZ(index(ii)));
            fprintf(wobsz,'\n%s %8.3f ',char(Atlas.id{index(ii)}),Atlas.DXYZ(index(ii)));

            fprintf(wpredx,'\n%s %8.3f ',char(Atlas.id{index(ii)}),Atlas.DXYZ(index(ii)));
            fprintf(wpredy,'\n%s %8.3f ',char(Atlas.id{index(ii)}),Atlas.DXYZ(index(ii)));
            fprintf(wpredz,'\n%s %8.3f ',char(Atlas.id{index(ii)}),Atlas.DXYZ(index(ii)));

            fprintf(wresx,'\n%s %8.3f ',char(Atlas.id{index(ii)}),Atlas.DXYZ(index(ii)));
            fprintf(wresy,'\n%s %8.3f ',char(Atlas.id{index(ii)}),Atlas.DXYZ(index(ii)));
            fprintf(wresz,'\n%s %8.3f ',char(Atlas.id{index(ii)}),Atlas.DXYZ(index(ii)));
            
        end
        

        fprintf(wobsx,'%12.8e ',dobs{1}(ii,jj));
        fprintf(wobsy,'%12.8e ',dobs{2}(ii,jj));
        fprintf(wobsz,'%12.8e ',dobs{3}(ii,jj));
        
        fprintf(wpredx,'%12.8e ',dpred{1}(ii,jj));
        fprintf(wpredy,'%12.8e ',dpred{2}(ii,jj));
        fprintf(wpredz,'%12.8e ',dpred{3}(ii,jj));
        
        fprintf(wresx,'%12.8e ',dres{1}(ii,jj));
        fprintf(wresy,'%12.8e ',dres{2}(ii,jj));
        fprintf(wresz,'%12.8e ',dres{3}(ii,jj));
        
        end
        
    end
    
end
    
fclose(rfid1);
fclose(rfid2);

if strcmp(argin,'write')==1
fclose(wobsx);fclose(wobsy);fclose(wobsz);
fclose(wpredx);fclose(wpredy);fclose(wpredz);
fclose(wresx);fclose(wresy);fclose(wresz);
end

%% Re-survey for plotting



% figure; plot((dpred{1}(:,1) - dobs{1}(:,1)) ./ dstd{1}(:,1),'*'); title('dBx - T1')
% figure; plot((dpred{1}(:,2) - dobs{1}(:,2)) ./ dstd{1}(:,2),'*'); title('dBx - T2')
% 
% figure; plot((dpred{2}(:,1) - dobs{2}(:,1)) ./ dstd{2}(:,1),'*'); title('dBy - T1')
% figure; plot((dpred{2}(:,2) - dobs{2}(:,2)) ./ dstd{2}(:,2),'*'); title('dBy - T2')

% Function to add random noise to ZTEM data

close all
clear all

workdir = 'C:\Projects\4017_ZTEM_lp_norm\Inversion\Dom';
cd(workdir);

fid=fopen('dpred_0.txt','rt');

out_name = '4017_base1_5freq_syn_5pc_noise_v2.obs';
writo = fopen(out_name,'a');

max_num_lines = 30000;
counter = 1;
flag = 0;
pct_noise = 0.1;
for ii=1:2
 line=fgets(fid); %gets next line
 
%     if ii==4
%        ndata = str2num(line)-1;
%         data = zeros(ndata,11); 
%     end
 fprintf(writo,'%s',line);     


end
fclose(writo);

line=fgets(fid);
counter = 1;

while line~=-1         	

    if isempty(findstr(line,'!frequency'))==0
        writo = fopen(out_name,'a');
        fprintf(writo,'%s',line);
        
        header1=fgets(fid);
        header2=fgets(fid);
        header3=fgets(fid);
        

        fclose(writo);
        
        counter = 0;
        data=zeros(1,11);
    end
    
    line=fgets(fid);
    
    while isempty(findstr(line,'!frequency'))==1 && length(line)>2
        counter = counter+1;
        data(counter,:)=str2num(line);
        line=fgets(fid);
        
    end
    
    
    if counter > 1 && (line(1))~=-1 
 %Add random error to each column 
        
        writo = fopen(out_name,'a');
        
               
        data(:,5) = pct_noise .* mean(data(:,4)) .* randn(counter,1);
        data(:,4) = data(:,4) + data(:,5);
        data(:,5) = abs(data(:,5));
        

        data(:,7) = pct_noise .* mean(data(:,6)) .* randn(counter,1);
        data(:,6) = data(:,6) + data(:,7);
        data(:,7) = abs(data(:,6))*0.02 + 0.0001;

        data(:,9) = pct_noise .* mean(data(:,8)) .* randn(counter,1);
        data(:,8) = data(:,8) + data(:,9);
        data(:,9) = abs(data(:,8))*0.02 + 0.0001;

        data(:,11) = pct_noise .* mean(data(:,10)) .* randn(counter,1);
        data(:,10) = data(:,10) + data(:,11);
        data(:,11) = abs(data(:,10))*0.02 + 0.0001;

        %% Keep only lines every 200 m
        
        fprintf(writo,'%s',num2str(sum(mod(data(:,2),200)==0)+1));
        fprintf(writo,'\n');
        fprintf(writo,'%s',header2);
        fprintf(writo,'%s',header3);
        fprintf(writo,'%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n',data((mod(data(:,2),200)==0),:)');
        fprintf(writo,'\n');
        fclose(writo);
        
        
%             while flag==0
%                 counter = counter +1;
%                 
%                 line=fgets(fid);
%                 
%                 if counter==2
%                 fprintf(writo,'%i\n',sum(mod(data(:,2),200)==0)+1);    
%                 elseif counter==1 || counter==3
%                 fprintf(writo,'%s',line);
%                 end
%                 
%                 
%                 
%                 if length(line) > 1
%                     
%                     if strcmp(line(1:2),'!X')==1                        
%                         data = zeros(ndata,11);
%                         flag = 1 ;
%                         counter = 1;
%                         fclose(writo);
%                     end
%                     
%                 else
%                     fprintf('Reach blank line. End of program\n')
%                     fclose(writo);
%                     break
%                     
%                 end
% 
%             end
            
            
    end
end
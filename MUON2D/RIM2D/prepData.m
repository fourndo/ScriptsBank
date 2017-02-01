 %getpath
%Compute N number of randomely generated ray paths.
%First determine if ray intersects a receiver
%If NO, ray path is only plotted.
%If YES, ray path populate a data matrix, then plotted

close all
clear all
addpath data/XH_RIM_755-759
addpath functions

DELIMITER = '\t';
HEADERLINES = 5;

% Import the file
in = importdata('Tx755_rx759_30m_down.xhs', DELIMITER, HEADERLINES);
RX=in.data(:,1)';
TX=in.data(:,2)';
Pha=in.data(:,3)';
Amp=in.data(:,4)';

marker=zeros(1,length(TX));


bin_range=2;


counter=1;

for ii=1:length(RX)
    
    if marker(ii)==0
    selector=( RX > (RX(ii)-bin_range)) .* (RX < ((RX(ii))+bin_range));
    
    Pha_bin(counter)=mean(Pha(selector==1));
    Amp_bin(counter)=mean(Amp(selector==1));
    RX_bin(counter)=(RX(ii));
    TX_bin(counter)=TX(ii);
    marker(selector==1)=1;
    
    counter=counter+1;
    end    
      
   
end

clear Rx Tx Pha Amp
% Import the file
in = importdata('Tx755_rx759_50m_down.xhs', DELIMITER, HEADERLINES);
RX=in.data(:,1)';
TX=in.data(:,2)';
Pha=in.data(:,3)';
Amp=in.data(:,4)';

marker=zeros(1,length(TX));

for ii=1:length(RX)
    
    if marker(ii)==0
    selector=( RX > (RX(ii)-bin_range)) .* (RX < ((RX(ii))+bin_range));
    
    Pha_bin(counter)=mean(Pha(selector==1));
    Amp_bin(counter)=mean(Amp(selector==1));
    RX_bin(counter)=(RX(ii));
    TX_bin(counter)=TX(ii);
    marker(selector==1)=1;
    
    counter=counter+1;
    end    
      
   
end

clear Rx Tx Pha Amp
% Import the file
in = importdata('Tx755_rx759_70m_down.xhs', DELIMITER, HEADERLINES);
RX=in.data(:,1)';
TX=in.data(:,2)';
Pha=in.data(:,3)';
Amp=in.data(:,4)';

marker=zeros(1,length(TX));

for ii=1:length(RX)
    
    if marker(ii)==0
    selector=( RX > (RX(ii)-bin_range)) .* (RX < ((RX(ii))+bin_range));
    
    Pha_bin(counter)=mean(Pha(selector==1));
    Amp_bin(counter)=mean(Amp(selector==1));
    RX_bin(counter)=(RX(ii));
    TX_bin(counter)=TX(ii);
    marker(selector==1)=1;
    
    counter=counter+1;
    end    
      
   
end

clear Rx Tx Pha Amp
% Import the file
in = importdata('Tx755_rx759_90m_down.xhs', DELIMITER, HEADERLINES);
RX=in.data(:,1)';
TX=in.data(:,2)';
Pha=in.data(:,3)';
Amp=in.data(:,4)';

marker=zeros(1,length(TX));

for ii=1:length(RX)
    
    if marker(ii)==0
    selector=( RX > (RX(ii)-bin_range)) .* (RX < ((RX(ii))+bin_range));
    
    Pha_bin(counter)=mean(Pha(selector==1));
    Amp_bin(counter)=mean(Amp(selector==1));
    RX_bin(counter)=(RX(ii));
    TX_bin(counter)=TX(ii);
    marker(selector==1)=1;
    
    counter=counter+1;
    end    
      
   
end

clear Rx Tx Pha Amp
% Import the file
in = importdata('Tx755_rx759_110m_down.xhs', DELIMITER, HEADERLINES);
RX=in.data(:,1)';
TX=in.data(:,2)';
Pha=in.data(:,3)';
Amp=in.data(:,4)';

marker=zeros(1,length(TX));

for ii=1:length(RX)
    
    if marker(ii)==0
    selector=( RX > (RX(ii)-bin_range)) .* (RX < ((RX(ii))+bin_range));
    
    Pha_bin(counter)=mean(Pha(selector==1));
    Amp_bin(counter)=mean(Amp(selector==1));
    RX_bin(counter)=(RX(ii));
    TX_bin(counter)=TX(ii);
    marker(selector==1)=1;
    
    counter=counter+1;
    end    
      
   
end

clear Rx Tx Pha Amp
% Import the file
in = importdata('Tx755_rx759_130m_down.xhs', DELIMITER, HEADERLINES);
RX=in.data(:,1)';
TX=in.data(:,2)';
Pha=in.data(:,3)';
Amp=in.data(:,4)';

marker=zeros(1,length(TX));

for ii=1:length(RX)
    
    if marker(ii)==0
    selector=( RX > (RX(ii)-bin_range)) .* (RX < ((RX(ii))+bin_range));
    
    Pha_bin(counter)=mean(Pha(selector==1));
    Amp_bin(counter)=mean(Amp(selector==1));
    RX_bin(counter)=(RX(ii));
    TX_bin(counter)=TX(ii);
    marker(selector==1)=1;
    
    counter=counter+1;
    end    
      
   
end

clear Rx Tx Pha Amp
% Import the file
in = importdata('Tx755_rx759_150m_down.xhs', DELIMITER, HEADERLINES);
RX=in.data(:,1)';
TX=in.data(:,2)';
Pha=in.data(:,3)';
Amp=in.data(:,4)';

marker=zeros(1,length(TX));

for ii=1:length(RX)
    
    if marker(ii)==0
    selector=( RX > (RX(ii)-bin_range)) .* (RX < ((RX(ii))+bin_range));
    
    Pha_bin(counter)=mean(Pha(selector==1));
    Amp_bin(counter)=mean(Amp(selector==1));
    RX_bin(counter)=(RX(ii));
    TX_bin(counter)=TX(ii);
    marker(selector==1)=1;
    
    counter=counter+1;
    end    
      
   
end

clear Rx Tx Pha Amp
% Import the file
in = importdata('Tx755_rx759_170m_down.xhs', DELIMITER, HEADERLINES);
RX=in.data(:,1)';
TX=in.data(:,2)';
Pha=in.data(:,3)';
Amp=in.data(:,4)';

marker=zeros(1,length(TX));

for ii=1:length(RX)
    
    if marker(ii)==0
    selector=( RX > (RX(ii)-bin_range)) .* (RX < ((RX(ii))+bin_range));
    
    Pha_bin(counter)=mean(Pha(selector==1));
    Amp_bin(counter)=mean(Amp(selector==1));
    RX_bin(counter)=(RX(ii));
    TX_bin(counter)=TX(ii);
    marker(selector==1)=1;
    
    counter=counter+1;
    end    
      
   
end

clear Rx Tx Pha Amp
% Import the file
in = importdata('Tx755_rx759_190m_down.xhs', DELIMITER, HEADERLINES);
RX=in.data(:,1)';
TX=in.data(:,2)';
Pha=in.data(:,3)';
Amp=in.data(:,4)';

marker=zeros(1,length(TX));

for ii=1:length(RX)
    
    if marker(ii)==0
    selector=( RX > (RX(ii)-bin_range)) .* (RX < ((RX(ii))+bin_range));
    
    Pha_bin(counter)=mean(Pha(selector==1));
    Amp_bin(counter)=mean(Amp(selector==1));
    RX_bin(counter)=(RX(ii));
    TX_bin(counter)=TX(ii);
    marker(selector==1)=1;
    
    counter=counter+1;
    end    
      
   
end

clear Rx Tx Pha Amp
% Import the file
in = importdata('Tx755_rx759_210m_down.xhs', DELIMITER, HEADERLINES);
RX=in.data(:,1)';
TX=in.data(:,2)';
Pha=in.data(:,3)';
Amp=in.data(:,4)';

marker=zeros(1,length(TX));

for ii=1:length(RX)
    
    if marker(ii)==0
    selector=( RX > (RX(ii)-bin_range)) .* (RX < ((RX(ii))+bin_range));
    
    Pha_bin(counter)=mean(Pha(selector==1));
    Amp_bin(counter)=mean(Amp(selector==1));
    RX_bin(counter)=(RX(ii));
    TX_bin(counter)=TX(ii);
    marker(selector==1)=1;
    
    counter=counter+1;
    end    
      
   
end
clear Rx Tx Pha Amp
% Import the file
in = importdata('Tx755_rx759_230m_down.xhs', DELIMITER, HEADERLINES);
RX=in.data(:,1)';
TX=in.data(:,2)';
Pha=in.data(:,3)';
Amp=in.data(:,4)';

marker=zeros(1,length(TX));

for ii=1:length(RX)
    
    if marker(ii)==0
    selector=( RX > (RX(ii)-bin_range)) .* (RX < ((RX(ii))+bin_range));
    
    Pha_bin(counter)=mean(Pha(selector==1));
    Amp_bin(counter)=mean(Amp(selector==1));
    RX_bin(counter)=(RX(ii));
    TX_bin(counter)=TX(ii);
    marker(selector==1)=1;
    
    counter=counter+1;
    end    
      
   
end
clear Rx Tx Pha Amp
% Import the file
in = importdata('Tx755_rx759_250m_down.xhs', DELIMITER, HEADERLINES);
RX=in.data(:,1)';
TX=in.data(:,2)';
Pha=in.data(:,3)';
Amp=in.data(:,4)';

marker=zeros(1,length(TX));

for ii=1:length(RX)
    
    if marker(ii)==0
    selector=( RX > (RX(ii)-bin_range)) .* (RX < ((RX(ii))+bin_range));
    
    Pha_bin(counter)=mean(Pha(selector==1));
    Amp_bin(counter)=mean(Amp(selector==1));
    RX_bin(counter)=(RX(ii));
    TX_bin(counter)=TX(ii);
    marker(selector==1)=1;
    
    counter=counter+1;
    end    
      
   
end
clear Rx Tx Pha Amp
% Import the file
in = importdata('Tx755_rx759_270m_down.xhs', DELIMITER, HEADERLINES);
RX=in.data(:,1)';
TX=in.data(:,2)';
Pha=in.data(:,3)';
Amp=in.data(:,4)';

marker=zeros(1,length(TX));

for ii=1:length(RX)
    
    if marker(ii)==0
    selector=( RX > (RX(ii)-bin_range)) .* (RX < ((RX(ii))+bin_range));
    
    Pha_bin(counter)=mean(Pha(selector==1));
    Amp_bin(counter)=mean(Amp(selector==1));
    RX_bin(counter)=(RX(ii));
    TX_bin(counter)=TX(ii);
    marker(selector==1)=1;
    
    counter=counter+1;
    end    
      
   
end


data=[RX_bin' TX_bin' Pha_bin' Amp_bin'];

save ('data/data','data');
save ('data/Tx','TX_bin');
save ('data/Rx','RX_bin');


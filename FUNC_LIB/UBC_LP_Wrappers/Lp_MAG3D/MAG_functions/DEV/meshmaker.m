function [dx,dy,dz]=meshmaker(meshfile)

header=6;

lineX=meshfile(1,1);
lineY=meshfile(2,1);
lineZ=meshfile(3,1);

%Space allocation for X, Y, Z vector position
dx=zeros(1,lineX);
dy=zeros(1,lineY);
dz=zeros(1,lineZ);

%Extract SouthWest top corner of grid

%Indices
index=header+1;
num=0;


while isnan(meshfile(index))==0 && index<=(ceil(lineX/3)*3+header)
    
    num=num+1;

    step=meshfile(index);
    dx(num)=step;
    index=index+1;
end

num=0;
index=ceil(lineX/3)*3+header+1;

while isnan(meshfile(index))==0 && index<=(ceil(lineX/3)*3+ceil(lineY/3)*3+header)
    
    num=num+1;
    step=meshfile(index);
    dy(num)=step;
    index=index+1;
end


num=0;
index=ceil(lineX/3)*3+ceil(lineY/3)*3+header+1;

while isnan(meshfile(index))==0 && index<=(ceil(lineX/3)*3+ceil(lineY/3)*3+ceil(lineZ/3)*3+header)
   
   num=num+1;
   step=meshfile(index);
   dz(num)=step;
   index=index+1;
end

%Create meshgrid
% [Xm,Ym,Zm]=meshgrid(Xstep,Ystep,Zstep);
% [dX,dZ,dY]=meshgrid(dx,dz,dy);
% 
% mesh_x=reshape(dX,lineX*lineY*lineZ,1);
% mesh_y=reshape(dY,lineX*lineY*lineZ,1);
% mesh_z=reshape(dZ,lineX*lineY*lineZ,1);

clear dX dZ dY
end
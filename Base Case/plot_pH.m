clear all
set(0, 'DefaultAxesFontSize', 18)


%Import data:
delimiterIn = ' ';
ZoneNum =7;
GridSize=3500;

if exist('co2d_conc.mat')==2
    load co2d_conc
else
    
    for i=1:ZoneNum
        %Reduce porosity
        headerlinesIn = 11+(i-1)*(GridSize+1);
        textdata = importdata('co2d_conc.dat',delimiterIn,headerlinesIn);
        zone=textdata.data;
        data(:,i)=zone(:,7);

        
    end
    
    save pH
end


%% Load data from the file TOUGHout_3D...

if exist('out.mat')==2 && exist('mesh_data.mat')==2
    load out
else
    [out times]=READ_DATA('flow.out');
    save out out times
end

%% load mesh data already saved by READ_DATA

load mesh_data

X=Coor{1,1};
Y=Coor{1,2};
Z=Coor{1,3};

pH=data(:,2);
x_mesh=out{1,1}(:,1); % also val{1,1} can be used from file mesh_data.mat
y_mesh=out{1,1}(:,2); % also val{1,2} can be used from file mesh_data.mat
z_mesh=out{1,1}(:,3); % also val{1,3} can be used from file mesh_data.mat


%% Plotting at fixed depth refining the mesh - PLOT XZ:
% we need to find in this block-by-block vector only the blocks that are
% at fixed depth. Let's say we meant to plot at z=-1810 (not that this
% number has to be one coordinate of the original TOUGH2 mesh)

xlin=linspace(min(X),max(X),1000); %linear spacing betwwen min(X) and max(X)
zlin=linspace(min(Z),max(Z),1000); %linear spacing between min(Y) and max(Y)

[Xcoord,Zcoord]=meshgrid(xlin,zlin);

var_image=griddata(x_mesh,z_mesh,pH,Xcoord,Zcoord);

figure
image(xlin,zlin,var_image,'Cdatamapping','scaled')
axis image
set(gca,'YDir','normal')
colorbar
xlabel('X (m)')
ylabel('Z (m)')
title('pH at 100 years')

print('pH_CO2_100.jpg','-djpeg','-r1200');



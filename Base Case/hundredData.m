clear all
set(0, 'DefaultAxesFontSize', 14)

%% Load data from the file flow.out

[out times]=READ_DATA('flow.out');
save out out times

%% load mesh data already saved by READ_DATA

load mesh_data
load volume

X=Coor{1,1};
Y=Coor{1,2};
Z=Coor{1,3};

x_mesh=out{1,1}(:,1); % also val{1,1} can be used from file mesh_data.mat
y_mesh=out{1,1}(:,2); % also val{1,2} can be used from file mesh_data.mat
z_mesh=out{1,1}(:,3); % also val{1,3} can be used from file mesh_data.mat

for i=1:6
    P=out{1,i}(:,4); % also val{1,3} can be used from file mesh_data.mat
    T=out{1,i}(:,5); % also val{1,3} can be used from file mesh_data.mat
    SG=out{1,i}(:,6); % also val{1,3} can be used from file mesh_data.mat
    XCO2=out{1,i}(:,10); % also val{1,3} can be used from file mesh_data.mat
    Pcap=out{1,i}(:,11); % also val{1,3} can be used from file mesh_data.mat
    Kred=out{1,i}(:,12); % also val{1,3} can be used from file mesh_data.mat
    DG=out{1,i}(:,13); % also val{1,3} can be used from file mesh_data.mat
    DL=out{1,i}(:,14); % also val{1,3} can be used from file mesh_data.mat
    
    %for j=1:1120
        %if XCO2(j)<0.001
            %XCO2(j)=0;
        %end
    %end
    
    %for j=1:1120
        %if XCO2(j)<0
            %XCO2(j)= -1 * XCO2(j);
        %end
    %end
    
    gas_co2(i)=vpa(sum( volume*0.25/1000000000.*SG.*DG,'all'),8)
    aqueous_co2(i)= vpa(sum( volume*0.25/1000000000.*(1-SG).*DL.*XCO2,'all'),8)
    
end


%% Plotting at fixed depth refining the mesh - PLOT XZ:
% we need to find in this block-by-block vector only the blocks that are
% at fixed depth. Let's say we meant to plot at z=-1810 (not that this
% number has to be one coordinate of the original TOUGH2 mesh)

xlin=linspace(min(X),max(X),1000); %linear spacing betwwen min(X) and max(X)
zlin=linspace(min(Z),max(Z),1000); %linear spacing between min(Y) and max(Y)
[Xcoord,Zcoord]=meshgrid(xlin,zlin);

figure
var_image=griddata(x_mesh,z_mesh,SG,Xcoord,Zcoord);
% h=figure('units','normalized','outerposition',[0 0 1 1]);
image(xlin,zlin,var_image,'Cdatamapping','scaled')
axis image
set(gca,'YDir','normal')
colorbar;
xlabel('Distance (m)')
ylabel('Depth (m)')
% set(gca,'FontSize',24)
% yticks([-24500 0])
% yticklabels({'-3500','0'})
title('Gas saturation at 500 years\newlineno minc, k_x=100 mD')
% yticks([0 200 400 600 800 1000])
% yticklabels({'-12','-10','-8','-6','-4','-2'})
hcb=colorbar;
% title(hcb,'Gas saturation')
print('SG500.jpg','-djpeg','-r1200');

figure
var_image=griddata(x_mesh,z_mesh,XCO2,Xcoord,Zcoord);
image(xlin,zlin,var_image,'Cdatamapping','scaled')
axis image
set(gca,'YDir','normal')
colorbar;
xlabel('Distance (m)')
ylabel('Depth (m)')
% yticks([-24500 0])
% yticklabels({'-3500','0'})
title('Dissolved CO2 at 500 years')
% yticks([0 200 400 600 800 1000])
% yticklabels({'-12','-10','-8','-6','-4','-2'})
hcb=colorbar;
% print('XCO2500.jpg','-djpeg','-r1200');

figure
var_image=griddata(x_mesh,z_mesh,DL,Xcoord,Zcoord);
image(xlin,zlin,var_image,'Cdatamapping','scaled')
axis image
set(gca,'YDir','normal')
colorbar;
xlabel('Distance (m)')
ylabel('Depth (m)')
% yticks([-24500 0])
% yticklabels({'-3500','0'})
title('Water density (kg/m^3) at 500 years')
% yticks([0 200 400 600 800 1000])
% yticklabels({'-12','-10','-8','-6','-4','-2'})
hcb=colorbar;
% print('DL500.jpg','-djpeg','-r1200');

figure
var_image=griddata(x_mesh,z_mesh,DG,Xcoord,Zcoord);
image(xlin,zlin,var_image,'Cdatamapping','scaled')
axis image
set(gca,'YDir','normal')
colorbar;
xlabel('Distance (m)')
ylabel('Depth (m)')
% yticks([-24500 0])
% yticklabels({'-3500','0'})
title('Gas density (kg/m^3) at 500 years')
% yticks([0 200 400 600 800 1000])
% yticklabels({'-12','-10','-8','-6','-4','-2'})
hcb=colorbar;
% print('DG500.jpg','-djpeg','-r1200');


figure
var_image=griddata(x_mesh,z_mesh,P,Xcoord,Zcoord);
% h=figure('units','normalized','outerposition',[0 0 1 1]);
image(xlin,zlin,var_image,'Cdatamapping','scaled')
axis image
set(gca,'YDir','normal')
colorbar;
xlabel('X (m)')
ylabel('Z (m)')
% set(gca,'FontSize',24)
% yticks([-24500 0])
% yticklabels({'-3500','0'})
title('Pressure (Pa) at 100 years')
% yticks([0 200 400 600 800 1000])
% yticklabels({'-12','-10','-8','-6','-4','-2'})
hcb=colorbar;
% title(hcb,'Gas saturation')
% print('SG500.jpg','-djpeg','-r1200');

figure
var_image=griddata(x_mesh,z_mesh,Pcap,Xcoord,Zcoord);
% h=figure('units','normalized','outerposition',[0 0 1 1]);
image(xlin,zlin,var_image,'Cdatamapping','scaled')
axis image
set(gca,'YDir','normal')
colorbar;
xlabel('X (m)')
ylabel('Z (m)')
% set(gca,'FontSize',24)
% yticks([-24500 0])
% yticklabels({'-3500','0'})
title('Capillary pressure (Pa) at 100 years')
% yticks([0 200 400 600 800 1000])
% yticklabels({'-12','-10','-8','-6','-4','-2'})
hcb=colorbar;
% title(hcb,'Gas saturation')
% print('SG500.jpg','-djpeg','-r1200');

figure
var_image=griddata(x_mesh,z_mesh,Kred,Xcoord,Zcoord);
% h=figure('units','normalized','outerposition',[0 0 1 1]);
image(xlin,zlin,var_image,'Cdatamapping','scaled')
axis image
set(gca,'YDir','normal')
colorbar;
xlabel('X (m)')
ylabel('Z (m)')
% set(gca,'FontSize',24)
% yticks([-24500 0])
% yticklabels({'-3500','0'})
title('Relative permeability at 100 years')
% yticks([0 200 400 600 800 1000])
% yticklabels({'-12','-10','-8','-6','-4','-2'})
hcb=colorbar;
% title(hcb,'Gas saturation')
% print('SG500.jpg','-djpeg','-r1200');
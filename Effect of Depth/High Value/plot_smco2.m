% clear all
set(0, 'DefaultAxesFontSize', 18)


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

x_mesh=out{1,1}(:,1); % also val{1,1} can be used from file mesh_data.mat
y_mesh=out{1,1}(:,2); % also val{1,2} can be used from file mesh_data.mat
z_mesh=out{1,1}(:,3); % also val{1,3} can be used from file mesh_data.mat


%Import data:
delimiterIn = ' ';
ZoneNum =8;
GridSize=1120;

for i=1:ZoneNum
    %Reduce porosity
    headerlinesIn = 11+(i-1)*(GridSize+1);
    textdata = importdata('co2d_min.tec',delimiterIn,headerlinesIn);
    zone=textdata.data;
    SMco2=zone(:,5);
    %for j=1:1120
        %if SMco2(j)<0.001
            %SMco2(j)=0;
        %end
    %end
    for j=1:1120
        if SMco2(j)<0
            SMco2(j)= -1 * SMco2(j);
        end
    end
    min_co2(i)=sum(SMco2.*volume/1000/1000000,'all');
    
end
%
%% Plotting at fixed depth refining the mesh - PLOT XZ:
% we need to find in this block-by-block vector only the blocks that are
% at fixed depth. Let's say we meant to plot at z=-1810 (not that this
% number has to be one coordinate of the original TOUGH2 mesh)

xlin=linspace(min(X),max(X),1000); %linear spacing betwwen min(X) and max(X)
zlin=linspace(min(Z),max(Z),1000); %linear spacing between min(Y) and max(Y)

[Xcoord,Zcoord]=meshgrid(xlin,zlin);

figure
var_image=griddata(x_mesh,z_mesh,SMco2,Xcoord,Zcoord);

image(xlin,zlin,var_image,'Cdatamapping','scaled')
axis image
set(gca,'YDir','normal')
hcb=colorbar
title('Mineralized CO_2 (kg/m^3 rock) at 100 years')

xlabel('X (m)')
ylabel('Depth (m)')
% yticks([-24500 0])
% yticklabels({'-3500','0'})
% print('SMCO2500.jpg','-djpeg','-r1200');

time=[0:10:50];

figure
box on
ax = gca;
ax.LineWidth = 1;
set(gca,'FontSize',12)
hold on
total_co2=gas_co2+aqueous_co2+min_co2(2:6);
plot(time,[0 total_co2],'-','LineWidth',1.5)
plot(time,[0 gas_co2],'--','LineWidth',1.5)
plot(time,[0 aqueous_co2],':','LineWidth',1.5)
plot(time,[0 min_co2(2:6)],'-.','LineWidth',1.5)
xlabel('Time (years)','FontWeight','bold')
ylabel('Amount of CO_2 (Mt)','FontWeight','bold')
ylim([0 4])
title('CO_2 injection, 25 ^oC')
legend('Total injected','Gaseous','Aqueous','Mineralized','location','northeast');
% print('TotalCO2_25C.jpg','-djpeg','-r1200');

save base

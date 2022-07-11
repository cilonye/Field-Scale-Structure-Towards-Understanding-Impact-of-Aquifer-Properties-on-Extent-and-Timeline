function [OUT,times]=READ_DATA(file,command)
% Read data from TOUGH2 output (MATLAB version for ext program)
% [OUT tempi]=READ_DATA(file,command,command2)-----------------------------
% This program doesn't read the secondary variables if included in the
% TOUGH2 output. It read just the normal output and eventually compute the
% fluxes in each grid block, since in the TOUGH2 output fluxes are given
% for the connection and not for the grid block.
% 
% INPUT variables
%    file -----> TOUGH2 output file name
%    command --> this command should be "Compute Flow" to read the flow
%                printout from the TOUGH output
%    command2 -> this command is only for Enhanced Coal Bed Methane version
%                of TOUGH2
% Variabili OUTPUT
%    OUT ------> Cell array 1xL where L is the number of the printout into 
%                the output. Each cell represent a printout in the TOUGH
%                output, and each cell is a NxM matrix where N is the
%                number of elements in the meshgrid and M is the number of
%                the variables such as pressure, temperature, etc..
%    times ----> Array 1xL of the printout time.
%--------------------------------------------------------------------------

t0=cputime;


if(nargin==3)
    if(strncmp(command,'Compute Flow',12) && strncmp(command2,'ECBM',4))
        compute_flux=1;
        option=1;
    else
        compute_flux=0;
        option=0;
        warning('Wrong command: the flow will not loaded')     
    end
elseif(nargin==2)
    if(strncmp(command,'Compute Flow',12))
        compute_flux=1;
        option=0;
    else
        compute_flux=0;
        option=0;
        warning('Wrong command: the flow will not loaded')
    end
elseif(nargin==1)
    compute_flux=0;
    option=0;
end

%%

% First of all this script read the MESH file using RMESH.m (see the script
% for more comments). 

if (compute_flux==1)
    if exist('mesh_conne.mat')==2
        load mesh_conne %The file MESH has been already read and all the 
                        %information, such as coordinates  
                        %and connection areas of the meshgrid, are stored 
                        %in the file "mesh_conne.mat"
        if(length(Coor{1,1})==1)
            disp(['Mesh ' num2str(length(Coor)-1) 'D'])
            disp(['with ' num2str(length(ID)) ' elements'])
            disp(['and ' num2str(length(conne)) ' connections'])
            disp('-------------------------------------------------------')
        else
            disp(['Mesh ' num2str(length(Coor)) 'D'])
            disp(['with ' num2str(length(ID)) ' elements'])
            disp(['and ' num2str(length(conne)) ' connections'])
            disp('-------------------------------------------------------')
        end
        x=val;
    else
        [Coor,ID,x,conne,ec,e,cc]=RMESH(); %Read the file mesh with all the 
                                           %possible options 
                                           %since the program READ_DATA has
                                           %been called with the option
                                           %"Compute_flow"
                                               
        if(length(Coor{1,1})==1)
            disp(['Mesh ' num2str(length(Coor)-1) 'D'])
            disp(['with ' num2str(length(ID)) ' elements'])
            disp(['and ' num2str(length(conne)) ' connections'])
            disp('-------------------------------------------------------')
        else
            disp(['Mesh ' num2str(length(Coor)) 'D'])
            disp(['with ' num2str(length(ID)) ' elements'])
            disp(['and ' num2str(length(conne)) ' connections'])
            disp('-------------------------------------------------------')
        end
        val=x;
        save mesh_conne Coor ID val conne ec e cc
    end
elseif (compute_flux==0)
    if exist('mesh_data.mat')==2
        load mesh_data %the file "mesh_data.mat" is a simplify version of 
                       %"mesh_conne.mat" since in this section the
                       %informations for flow calculation are useless
                       
        if (length(Coor{1,1})==1)
            disp(['Mesh ' num2str(length(Coor)-1) 'D'])
            disp(['with ' num2str(length(ID)) ' elements'])
            disp('-------------------------------------------------------')
        else
            disp(['Mesh ' num2str(length(Coor)) 'D'])
            disp(['with ' num2str(length(ID)) ' elements'])
            disp('-------------------------------------------------------')
        end
        x=val;
    else
        [Coor,ID,x]=RMESH(); %Read only the ELEME part of the MESH file. 

        if (length(Coor{1,1})==1)
            disp(['Mesh ' num2str(length(Coor)-1) 'D'])
            disp(['with ' num2str(length(ID)) ' elements'])
            disp('-------------------------------------------------------')
        else
            disp(['Mesh ' num2str(length(Coor)) 'D'])
            disp(['with ' num2str(length(ID)) ' elements'])
            disp('-------------------------------------------------------')
        end
        val=x;
        save mesh_data Coor ID val
    end
end
if(length(x{1,2})==length(x{1,1}))
    nxyz=3;
    coor_mesh=x;
elseif(length(Coor{1,1})==1)
    nxyz=1;
    coor_mesh{1,1}=x{1,1};
    coor_mesh{1,2}=x{1,3};
else
    nxyz=2;
    coor_mesh{1,1}=x{1,1};
    coor_mesh{1,2}=x{1,3};
end

%%

disp(['Loading data from "' file '".........'])
disp('')

if(~(exist(file)==2))
    error(['The file: "' file '" does not exist'])
end

iprnt = 0;
fid=fopen(file);

while(iprnt>=0)

    i11 = 1;
    lc(1:140)=' ';
    line(1:140)=' ';
    
    %Read through the output file and find the string "TOTAL TIME". The
    %program finish when the string "END OF TOUGH2 SIMULATION RUN" is
    %found.

    while (~(strncmp(lc(i11:i11+9),'total time',10)))
        line=fgets(fid);
        line(1)=' ';
        [lc,i11,l,i1l]=INSPECT(line,lc); % INSPECT function evaluate the 
                                         % position "i11" of the first not null 
                                         % character in "line", the
                                         % position of the last non-space 
                                         % character "l" and the position
                                         % "i1l" of the last character of
                                         % the first word. "lc" is the line
                                         % in lower-case
                                         
        if(length(line)>40 & (strncmp(line(2:29),'END OF TOUGH2 SIMULATION RUN',28) || strncmp(line(2:25),'END OF TOUGH2 SIMULATION',24)))
            fclose(fid);
            t=cputime-t0;
            disp('')
            disp('-------------------------------------------------------')
            disp('')
            disp(['Loading time: ' num2str(t)])
%             if (nxyz==1)
%                 for i=1:iprnt
%                     if (compute_flux==0)
%                         OUT{1,i}=[coor_mesh{1,2} v{1,i}];
%                     else
%                         OUT{1,i}=[coor_mesh{1,2} v{1,i} v2{1,i}];
%                     end
%                 end
%             elseif (nxyz==2)
%                 for i=1:iprnt
%                     if (compute_flux==0)
%                         OUT{1,i}=[coor_mesh{1,1:2} v{1,i}];
%                     else
%                         OUT{1,i}=[coor_mesh{1,1:2} v{1,i} v2{1,i}];
%                     end
%                 end
%             elseif(nxyz==3)
%                 for i=1:iprnt
%                     if (compute_flux==0)
%                         OUT{1,i}=[coor_mesh{1,1:3} v{1,i}];
%                     else
%                         OUT{1,i}=[coor_mesh{1,1:3} v{1,i} v2{1,i}];
%                     end
%                 end
%             end
            return
        end
    end

    line=fgets(fid);
    iprnt = iprnt + 1;

    disp(['Found printout(' num2str(iprnt) '):' line(3:20)])
    times(iprnt)=str2double(line(3:13));

    i11 = 1;

    while(~(strncmp(lc(i11:i11+3),'elem',4)))
        line=fgets(fid);
        line(1:1)=' ';
        [lc,i11,l,i1l]=INSPECT(line,lc);
    end

    i2l=i1l;

    [lc(i2l+1:l),ii1,l9,iil]=INSPECT(lc(i2l+1:l),lc(i2l+1:l));
    ii1=ii1+i2l;
    iil=iil+i2l;

    il=iil;
    nn=0;
    nd=0;

% determine which columns contain the numbers

    while(il<l)
        [lc(il+1:l),i1,l9,ll]=INSPECT(lc(il+1:l),lc(il+1:l));
        i1=i1+il;
        il=ll+i1;
        nn=nn+1;
        c1(nn)=i1;
        cl(nn)=il;
        uv(nn)=0;
        nd=nd+1;
    end

    l0=0;
    nu2=-1;
    index=1;

    var_temp=fgets(fid);

    if(length(var_temp)>5 && (var_temp(18)=='(' || var_temp(18)==')'))
        fgets(fid);
    end
    
    clear var_temp
    
    v=zeros(length(ID),nn);
         
%     a line of @'s ends this section
%     blank lines are ignored as are the page break headings
    

    while(~(strncmp(line(2:21),'@@@@@@@@@@@@@@@@@@@@',20)))
        line=fgets(fid);
        line(1:1)=' ';
        if (length(line)<5)
            for i=1:2
                temp=fgets(fid);
                if(length(temp)>20 & strncmp(temp(2:21),'@@@@@@@@@@@@@@@@@@@@',20))
                    line_temp=temp;
                    clear temp
                end
                clear temp
            end
            if(exist('line_temp'))
                line=line_temp;
                clear line_temp
            else
                line=fgets(fid);
                if (length(line)>5 && line(18)==')') %modified for ECO2N Berkeley 03/31/11
                    line=fgets(fid);
                    if length(line)<5
                        line=fgets(fid);
                    end
                elseif(length(line)<5)
                    line=fgets(fid);
                end
            end
        end

        if(l0==0)
            [lc,l0,l,k]=INSPECT(line(ii1:end),lc);
            l=l+ii1-1;
            if(cl(nn)<l)
                cl(nn)=l;
            end
            cl(nn+1)=l+1;
            iil=ii1+k-1;
            il=iil;
            for i=1:nn
                i1=il;
                il=cl(i+1);
                if(il>length(line))
                    [lc,l0,l,k]=INSPECT(line(i1+1:end),lc);
                else
                    [lc,l0,l,k]=INSPECT(line(i1+1:il),lc);
                end
                il=i1+k;
                if(k>40)
                    k=40;
                end
                i1=il-k;
                c1(i)=i1+1;
                cl(i)=il;
            end
        end
        
        % "nn" represents the number of parameters found in the printout,
        % such as pressure, temperature, etc..., generally nn=10. Variables
        % "c1(j)" and "cl(j)" represent the position of the first and last
        % character respectively for the parameter "j" in the read line

        if(~(strncmp(line(2:21),'@@@@@@@@@@@@@@@@@@@@',20)))
            for j=1:nn
                v(index,j)=str2double(line(c1(j):cl(j)));
                %% Berkeley 02/08/10
                % when the esponent of some values is very small, such as
                % smaller than -100, the value in the TOUGH2 output will be
                % something like 3.4567-100 instead of 3.4567E-100. For
                % this very small number a conversion in MATLAB with
                % str2double will produce a NaN. Then the values will be
                % set to 0. Comment the following 3 lines if you have
                % already some NaN in the TOUGH2 output
                if isnan(v(index,j))
                    v(index,j)=0;
                end
                %%
            end
            index=index+1;
        end
    end

% This part of the program is for the calculation of the fluxes in each
% grid block

    if(compute_flux==1)
 
        i11 = 1;
        lc(1:140)=' ';
        
        %find the line where the flow printout begin

        while(~(strncmp(lc(i11:i11+3),'elem',4)))
            line=fgets(fid);
            line(1:1)=' ';
            [lc,i11,l,i1l]=INSPECT(line,lc);
        end
        
        %this part is for ECBM version of TOUGH2 only
        if(option==1)
            line=fgets(fid);
            line(1:1)=' ';
            [lc,i11,l,i1l]=INSPECT(line,lc);
            while(~(strncmp(lc(i11:i11+4),'elem1',5)))
                line=fgets(fid);
                line(1:1)=' ';
                [lc,i11,l,i1l]=INSPECT(line,lc);
            end
        end
        %

        i2l=i1l;
        [lc(i2l+1:l),ii1,l9,iil]=INSPECT(lc(i2l+1:l),lc(i2l+1:l));
        ii1=ii1+i2l;
        iil=iil+i2l;
        while(~(strncmp(lc(ii1:iil),'index',5)))
            i21=ii1;
            i2l=iil;
            if(strncmp(lc(ii1:iil),'elem2',5))
                [lc(i2l+1:l),ii1,l9,iil]=INSPECT(lc(i2l+1:l),lc(i2l+1:l));
                ii1=ii1+i2l;
                iil=iil+i2l;
            end
        end
       

        il=iil;
        nn=0;
        nd=0;

        clear c1 cl uv
        
        % As in the first part of the program "nn", "c1", "cl" are useful
        % to get the number of fluxes, such as heat, gas flow, liquid flow,
        % and velocity, and their position in "line".

        while(il<l)
            [lc(il+1:l),i1,l9,ll]=INSPECT(lc(il+1:l),lc(il+1:l)); %"lc" is the "line" in lower-case
            i1=i1+il;
            il=ll+i1;
            nn=nn+1;
            c1(nn)=i1;
            cl(nn)=il;
            %uv(nn)=0;
            %nd=nd+1;
            % Added Berkeley 01/26/10 to solve problem (maybe)
            uv(nn)=1; % "uv" is a control variable, used to determine if we
                      % need to divide for the connection surface or not!
            if(strncmp(lc(i1:il+2),'vel',3))
                uv(nn)=-1;
            end
%            l9=il+2-i1;
            % v1(lv+1:lv+9)=line(i1-1:il); %header: I don't need in matlab!
%            i2=lv+2;
            % header: not for matlab!
%             if(itype==5 && strncmp(v1(i2+3:i2+3),'(',1))
%                 v1(i2+1:lv+l9)=v1(i2+4:lv+l9);
%                 l9=l9-3;
%                 if(l9>6 && strncmp(v1(i2+4:i2+4),')',1))
%                     v1(i2+2:lv+9)=v1(i2+5:lv+l9);
%                     l9=l9-3;
%                 end
%             elseif(itype==6 && strncmp(v1(i2+1:i2+1),'(',1))
%                 v1(i2+1:lv+l9)=v1(i2+2:lv+l9);
%                 l9=l9-2;
%             end
            %....other stuff for header (I hope!)            
        end
%% Added Berkeley 01/26/10
%         nd=nn-nd;
%         if(itype==5)
%             nd=nd*nxyz;
%         end
%         for i=1:nd
%             for nu1=1:ne
%                 v(i+nv,nu1)=0;
%             end
%         end
%         i=nv;
%         nv=nv+nd;
%         nd=i;

%%
        
        l0=0;
        nu2=-1;
        index=1;

        fgets(fid);

        if(nxyz==1)
            lenV2=nn;
        elseif(nxyz==2)
            lenV2=nn*2;
        elseif(nxyz==3)
            lenV2=nn*3;
        end

        v2=zeros(length(ID),lenV2);
        
        inde=0;
        jjjj=0;
%        for I=1:55
        while(~(strncmp(line(2:21),'@@@@@@@@@@@@@@@@@@@@',20)))
%            fgets(fid);
            line=fgets(fid);
            jjjj=jjjj+1;
             if (length(line)>=27 & strncmp(line(29),'(',1))
                 fgets(fid);
                 line=fgets(fid);
             end
             
             if (length(line)>=24 & strncmp(line(23),'(',1))
                 fgets(fid);
                 line=fgets(fid);
             end
            line(1:1)=' ';
            if (length(line)<5)
%                 for i=1:2
%                     temp=fgets(fid);
%                     if(length(temp)>20 & strncmp(temp(2:21),'@@@@@@@@@@@@@@@@@@@@',20))
%                         line_temp=temp;
%                         clear temp
%                     end
%                     clear temp
%                 end
                if(exist('line_temp'))
                    line=line_temp;
                    clear line_temp
                else
                    line=fgets(fid);
                    while(length(line)<5 | ((length(line)>=27 & strncmp(line(29),'(',1))) | ((length(line)>=24 & strncmp(line(23),'(',1))))
                        line=fgets(fid);
                    end
                end
            end
            
            % fluxes in the TOUGH2 output are given for the connection of
            % two elements in the gridmesh. "line(i11:i1l)" and
            % "line(i2i:i2l)" represent the identifier of such elements,
            % and "nu1" and "nu2" are two numbers indentifying the position
            % in the meshgrid, i.e. in the ID array producted with RMESH.
            
            
%            disp(line(i11:i1l))
            
            nu1=strmatch(line(i11:i1l),ID,'exact');
            nu2=strmatch(line(i21:i2l),ID,'exact');
            
            
                
            if ((isempty(nu1) | isempty(nu2) | l0==0) & strncmp(line(i11:i1l-1),'  ',2))
                nu1=strmatch(line(i11+2:i1l+2),ID,'exact');
                nu2=strmatch(line(i21+2:i2l+2),ID,'exact');
                i11=i11+2;
                i1l=i1l+2;
                i21=i21+2;
                i2l=i2l+2;
            end
            
            
            %% Berkeley 8/10/2011
            if ((isempty(nu1) | isempty(nu2) | l0==0) & ~strncmp(line(i11-1:i11),' ',1))
                nu1=strmatch(line(i11-1:i1l-1),ID,'exact');
                nu2=strmatch(line(i21-1:i2l-1),ID,'exact');
                i11=i11-1;
                i1l=i1l-1;
                i21=i21-1;
                i2l=i2l-1;
            end
            
            %% Berkeley 3/20/2012
            if ((isempty(nu1) | isempty(nu2) | l0==0) & ~strncmp(line(i11-1:i11),'  ',2))
                nu1=strmatch(line(i11-1:i1l-1),ID,'exact');
                nu2=strmatch(line(i21-1:i2l-1),ID,'exact');
                i11=i11-1;
                i1l=i1l-1;
                i21=i21-1;
                i2l=i2l-1;
            end
            
            
            % here the distance between the 2 elements wil be computed, in
            % the x, y, and z direction respectively. The fluxes through a
            % connection will be divided by these distances to have the
            % flux in a given element
            

            if(nu1~=0 & nu2~=0)
                dx=x{1,1}(nu1)-x{1,1}(nu2);
                dz=x{1,3}(nu1)-x{1,3}(nu2);
                if (nxyz==1)
                    denom1=2*dz;
                    dz=dz/denom1;
                elseif (nxyz==3)
                    dy=x{1,2}(nu1)-x{1,2}(nu2);
                    denom1=2*sqrt(dx*dx+dy*dy+dz*dz);
                    dx=dx/denom1;
                    dy=dy/denom1;
                    dz=dz/denom1;
                else
                    denom1=2*sqrt(dx*dx+dz*dz);
                    dx=dx/denom1;
                    dz=dz/denom1;
                end
                
                % connection has a identifying number. This number
                % is in the arrays "ec", "e", and "cc" produced by RMESH.
                
                nc=ec(nu1); 
                while(nc~=0 & nu2~=e(2,nc))
                    i=(nu1-e(1,nc))/(e(2,nc)-e(1,nc));
                    nc=cc(i+1,nc);
                end
                if (nc==0)
                    disp('no connection')
                    return
                end
             
                     
                %denom=conne(nc);
                %denom=1.0;
        
                % Here the position of the numbers for the fluxes will be
                % computed. As previously "c1(j)" and "cl(j)" represent the
                % position of th first and the last character in the line
                % for the "j"-flux.


                if(l0==0)
                    [lc,l0,l,k]=INSPECT(line(ii1:end),lc);
                    l=l+ii1-1;
                    if(cl(nn)<l)
                        cl(nn)=l;
                    end
                    cl(nn+1)=l+1;
                    iil=ii1+k-1;
                    il=iil;
                    for i=1:nn
                        i1=il;
                        il=cl(i+1);
                        if(il>length(line))
                            [lc,l0,l,k]=INSPECT(line(i1+1:end),lc);
                        else
                            [lc,l0,l,k]=INSPECT(line(i1+1:il),lc);
                        end
                        il=i1+k;
                        if(k>40)
                            k=40;
                        end
                        i1=il-k;
                        c1(i)=i1+1;
                        cl(i)=il;
                    end
                end

                if(~(strncmp(line(2:21),'@@@@@@@@@@@@@@@@@@@@',20)))
                    k=1;
                    for j=1:nn
                        %% Berkeley 01/26/10
                        denom2=1;
                        % if the fluxes aren't velocity we need to divide
                        % them by the connection area. Connection areas are
                        % in the array "conne" produced by RMESH.
                        if uv(j)>0
                            denom2=conne(nc);
                        end
                        %%
                        pv=(str2double(line(c1(j):cl(j))))/denom2;
                        %% Berkeley 02/08/10
                        % for very small number with exponent smaller than
                        % e-100....comment the next 3 lines if you have NaN
                        % in the fluxes printout of the original TOUGH2
                        % output
                        
                        if isnan(pv)
                            pv=0;
                        end
                        %%
                          
                        if(nxyz==1)
                            v2(nu1,k)=v2(nu1,k)+dz*pv;
                            v2(nu2,k)=v2(nu2,k)+dz*pv;
                            k=k+1;
                        elseif(nxyz==2)
                            v2(nu1,k)=v2(nu1,k)+dx*pv;
                            v2(nu2,k)=v2(nu2,k)+dx*pv;
                            v2(nu1,k+1)=v2(nu1,k+1)+dz*pv;
                            v2(nu2,k+1)=v2(nu2,k+1)+dz*pv;
                            k=k+2;
                        elseif(nxyz==3)
                            v2(nu1,k)=v2(nu1,k)+dx*pv;
                            v2(nu2,k)=v2(nu2,k)+dx*pv;
                            v2(nu1,k+1)=v2(nu1,k+1)+dy*pv;
                            v2(nu2,k+1)=v2(nu2,k+1)+dy*pv;
                            v2(nu1,k+2)=v2(nu1,k+2)+dz*pv;
                            v2(nu2,k+2)=v2(nu2,k+2)+dz*pv;
                            k=k+3;
                        end
                    end
                end
            end
        end
%         disp(['I= ' num2str(I)])
%         disp(line)
%         disp(['n1= ' num2str(nu1)])
%         disp(['n2= ' num2str(nu2)])
%         disp(dx)
%         disp(['dy= ' num2str(dy)])
%         disp(x{1,2}(nu1))
%         disp(x{1,2}(nu2))
%         disp(x{1,3}(nu1))
%         disp(x{1,3}(nu2))
%         disp(['dz= ' num2str(dz)])
%        end
    end
    
    % Finally the OUT is here assembled....the final output is a cell array
    % 1xL where L is the number of printout in the TOUGH2 output. Each cell
    % of OUT is a matrix with NxM elements, where N is the number of
    % elements in the meshgrid and M the number of the parameters. M
    % depends on the meshgrid dimension (1D,2D,or 3D) and on the flow
    % printout if wanted. Flow printout (v2) dimensions depend also on the
    % dimension. If "nn" is the number of fluxes in the TOUGH2 output, then
    % length(v2) will be "nn" for 1D, "2 x nn" for 2D, "3 x nn" for 3D
    
    if (nxyz==1)  
        if (compute_flux==0)
            OUT{1,iprnt}=[coor_mesh{1,2} v];
        else
            OUT{1,iprnt}=[coor_mesh{1,2} v v2];
        end
    elseif (nxyz==2)
        if (compute_flux==0)
            OUT{1,iprnt}=[coor_mesh{1,1:2} v];
        else
            OUT{1,iprnt}=[coor_mesh{1,1:2} v v2];
        end
    elseif(nxyz==3)
        if (compute_flux==0)
            OUT{1,iprnt}=[coor_mesh{1,1:3} v];
        else
            OUT{1,iprnt}=[coor_mesh{1,1:3} v v2];
        end
    end
end

fclose(fid);

t=t0-cputime
end

function [lc,l0,l,k0]=INSPECT(line,lc)
% Find in a line the first non-space character and the first word length 
% [lc,l0,l,k0]=INSPECT(line)-----------------------------------------------
% 
% INPUT Variables
%   line--> line (string) to analyze
%   lc----> empty string
% OUTPUT Variables
%   lc----> variable "line" in lower-case
%   l0----> position in line of the first non-space character
%   l-----> position in line of the last non-space character
%   k0----> position in line of the last character of the first word
% -------------------------------------------------------------------------

len1=length(line);
len2=length(lc);

if(len1>len2) 
    len1=len2;
end

%l0=0;
k0=0;
%l=1;

lc=lower(line);

vec=not(logical(isspace(line)));
vec(length(vec)+1)=0;

l0=find(vec,1,'first');
l=find(vec,1,'last');

i=l0;

while (vec(i)==1)
    k0=i;
    i=i+1;    
end


% for i=1:len1
%     c=line(i);
%     if(c~=' ')
% %        ic=double(c)+32;
% %        if(c>='A' && c<='Z')
% %            c=char(ic);
% %        end
%         if(l0==0)
%             l0=i;
%         end
%         if(k0==0)
%             k0=-1;
%         end
%         l=i;
%     elseif(k0<0)
%         k0=l;
%     end
% %    lc(i)=c;
% end



if (len2>len1)
    for i=len1+1:len2
        lc(i)=' ';
    end
end
if(l0==0)
    l0=1;
end
if(k0<1)
    k0=l;
end
end


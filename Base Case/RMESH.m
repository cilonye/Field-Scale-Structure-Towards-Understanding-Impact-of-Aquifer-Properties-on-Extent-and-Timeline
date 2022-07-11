function [coor,ID,coor_mesh,conne,ec,e,cc]=RMESH()
% File MESH reader
% [coor,ID,coor_mesh,conne]=RMESH()----------------------------------------
%              
% OUTPUT Variables
% coor -------> provides the meshgrid coordinates. It is a 2 cell array for
%               2D mesh, and a 3 cell array for 3D mesh. A cell array 
%               provides a storage mechanism for dissimilar kinds of data.
% ID ---------> Array of the elements name (5 char), with a length equal to
%               the mesh length.
% coor_mesh --> It is a 3 cell array, providing the coordinate X, Y, Z for
%               each meshgrid element. Each matrix of the cell array is
%               1xlength(ID).
% conne ------> Array containing the connection areas. It is useful to
%               compute the fluxes.
% ec,e,cc-----> hash table
% -------------------------------------------------------------------------

t0=cputime;

%% Analysis of the grid mesh

disp('-------------------------------------------------------')
disp('Reading file MESH.....')
disp('')

fid=fopen('MESH.txt');

if(fid<0)
    error('The file: "MESH" does not exist')
end

C{1,1}=fgets(fid);

if ((strncmp(C{1,1}(1,1:5),'ELEME',5)) || (strncmp(C{1,1}(1,1:5),'eleme',5)))
else
    error('no eleme in MESH')
end

%C1=textscan(fid,'%5c%61c',1);

%% Berkeley 03/23 for In Salah

C1{1,1}='     ';
i=0;

while ~((strncmp(C1{1,1}(1,1:5),'CONNE',5)) || (strncmp(C1{1,1}(1,1:5),'conne',5)))
    temp=fgets(fid);
    if length(temp)<5
        temp(length(temp):5)=' ';
    end
    C1{1,1}(1,1:5)=temp(1:5);
    i=i+1;
end

lmesh=i-2;

%%


% for i=1:length(C1{1,1})
%     if ((strncmp(C1{1,1}(i,1:5),'CONNE',5)) || (strncmp(C1{1,1}(i,1:5),'conne',5)))
%         lmesh=i-1;
%     end
% end

fclose(fid);
clear C C1;
fid=fopen('MESH.txt');

fgets(fid);
for i=1:lmesh
    temp=fgets(fid);
    C1{1,1}(i,1:5)=temp(1:5);
    C1{1,2}(i,1:45)=temp(6:50);
    C1{1,3}(i,:)=temp(51:80);
end

for i=1:lmesh
    ID(i,:)=C1{1,1}(i,:);
    xyz(i,:)=C1{1,3}(i,:);
end

if (strncmp(xyz(1,11:20),'          ',10) || (strncmp(xyz(1,11:20),'0.0000E+00',10) && strncmp(xyz(2,11:20),'0.0000E+00',10)));
    x(:,1)=str2num(xyz(:,1:10));
    y=[];
    z(:,1)=str2num(xyz(:,21:30));
    k=1;
    j=1;
    x2=sort(x);
    X(1)=x2(1);
    z2=sort(z);
    Z(1)=z2(1);
    for i=2:lmesh
        if (z2(i)~=z2(i-1))
            k=k+1;
            Z(k)=z2(i);
        end
        if (x2(i)~=x2(i-1))
            j=j+1;
            X(j)=x2(i);
        end
    end
    
    if Z<0
        Z=-Z;
        Z=sort(Z);
        Z=-Z;
    end
    if X<0
        X=-X;
        X=sort(X);
        X=-X;
    end
    
    coor{1,1}=X;
    coor{1,2}=Z;
        
else
    x(:,1)=str2num(xyz(:,1:10));
    y(:,1)=str2num(xyz(:,11:20));
    z(:,1)=str2num(xyz(:,21:30));
    k=1;
    j=1;
    n=1;
    x2=sort(x);
    X(1)=x2(1);
    y2=sort(y);
    Y(1)=y2(1);
    z2=sort(z);
    Z(1)=z2(1);
    for i=2:lmesh    
        if (z2(i)~=z2(i-1))
            k=k+1;
            Z(k)=z2(i);
        end
        if (x2(i)~=x2(i-1))
            j=j+1;
            X(j)=x2(i);
        end
        if (y2(i)~=y2(i-1))
            n=n+1;
            Y(n)=y2(i);
        end
    end
    
    if Z<0
        Z=-Z;
        Z=sort(Z);
        Z=-Z;
    end
    if Y<0
        Y=-Y;
        Y=sort(Y);
        Y=-Y;
    end
    if X<0
        X=-X;
        X=sort(X);
        X=-X;
    end
    
    coor{1,1}=X;
    coor{1,2}=Y;
    coor{1,3}=Z;
end

coor_mesh{1,1}=x;
coor_mesh{1,2}=y;
coor_mesh{1,3}=z;

[hsh,ec]=hash(ID); %hash(ID) is a function helpful to index the elements


%% Analysis of the connection!

if (nargout>3)
    fgets(fid);
    wrd=fgets(fid);

    if ((strncmp(wrd(1:5),'CONNE',5)) || (strncmp(wrd(1:5),'conne',5)))
    else
        error('no conne in MESH')
    end

    num_conne=0;

    line=fgets(fid);

    while(length(line)>9)
        num_conne=num_conne+1;
        a(num_conne)=str2double(line(51:60)); %Connection Area
        if (nargout>4)
            wrd=line(1:5);
            nu1=strmatch(wrd,ID,'exact');
%            nu1=ihash(wrd,ID);
            wrd2=line(6:10);
            nu2=strmatch(wrd2,ID,'exact');
%            nu2=ihash(wrd2,ID);
            e(1,num_conne)=nu1;
            e(2,num_conne)=nu2;
            cc(1,num_conne)=ec(nu1);
            cc(2,num_conne)=ec(nu2);
            ec(nu1)=num_conne;
            ec(nu2)=num_conne;
        end
        line=fgets(fid);
    end
    conne=a;
end

fclose(fid);

t=cputime-t0;

disp('')
disp(['Reading time: ' num2str(t)])
disp('-------------------------------------------------------')

end
%% FUNCTIONS

function [hsh,h,ind]=hash(ID)

lmesh=length(ID);
ind(1:lmesh)=1:1:lmesh;
h=zeros(lmesh,1);

for j=1:lmesh
    if (ID(j,4)=='0') 
        ID(j,4)=' ';
    end

    w1=ID(j,:);
    n=0;
    for i=1:5
        n=n+i*double(w1(i));
    end
    h(j)=double(int16(n/256))*256 + 1;
end

hsh=zeros(int32(max(h(j))+1),1);

i1=0;
for i=1:(max(h(j))+1)
    i2=i1+1;
    for j=i2:lmesh
        k=ind(j);
        if (h(k)==i)
            h(k)=0;
            i1=i1+1;
            ind(j)=ind(i1);
            ind(i1)=k;
        end
    end
    hsh(i)=i1;
end
end

function numb=ihash(wrd,ID)

ne=length(ID);

w1=wrd;
if (w1(4)=='0')
    w1(4)==' '
end

for i=1:ne
    if(strncmp(ID(i,:),w1,5))
        numb=i;
        return
    end
end
numb=0;

end

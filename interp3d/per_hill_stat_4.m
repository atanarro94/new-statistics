close all; 
clear all; 
clc
double precision;
format long;

%Load file containing number of points per core
load 'history2.txt'
ncores=history2(end,1)+1; % number of cores

%Read interpolating mesh
fname_x           = 'ZSTAT/x.fort';
[fid_x,message_x] = fopen(fname_x,'r','ieee-le');
hdr_x             = fread(fid_x,1,'int32');
dum1              = fread(fid_x,3,'int32');
npoints=dum1(1);
x            = fread(fid_x,npoints,'*float64');
fclose(fid_x);

fname_y           = 'ZSTAT/y.fort';
[fid_y,message_y] = fopen(fname_y,'r','ieee-le');
hdr_y             = fread(fid_y,1,'int32');
dum1              = fread(fid_y,3,'int32');
npoints=dum1(1);
y            = fread(fid_y,npoints,'*float64');
fclose(fid_y);

fname_z           = 'ZSTAT/z.fort';
[fid_z,message_z] = fopen(fname_z,'r','ieee-le');
hdr_z             = fread(fid_z,1,'int32');
dum1              = fread(fid_z,3,'int32');
npoints=dum1(1);
z            = fread(fid_z,npoints,'*float64');
fclose(fid_z);

% figure(1)
% plot3(z,x,y,'ob'); hold on;
% axis equal
% xlabel('z')
% ylabel('x')
% zlabel('y')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LOAD INTERPOLATED DATA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for field_number=1:22
%U
%%%%%
fname=['ZSTAT/U',num2str(field_number,'%2.2d')];

[fid,message] = fopen(fname,'r','ieee-le');
hdr           = fread(fid,1,'int32')      ;
F             = fread(fid,hdr,'*char')'   ;
dum5          = fread(fid,1,'*float64')   ;
time          = fread(fid,1,'*float64')   ;
sfn           = fread(fid,1,'int32')      ;
dum6          = fread(fid,1,'*float64')   ;

aa=fread(fid,Inf,'*float64'); 

%Assign points taking into account gaps among processors
for i=1:ncores
    bb(1+sum(history2(1:i-1,2)):sum(history2(1:i,2)),1)=aa(1+sum(history2(1:i-1,2))+(i-1):sum(history2(1:i,2))+(i-1));
end

UD{field_number}=bb; 

fclose(fid);

%W
%%%%%
fname=['ZSTAT/W',num2str(field_number,'%2.2d')];

[fid,message] = fopen(fname,'r','ieee-le');
hdr           = fread(fid,1,'int32')      ;
F             = fread(fid,hdr,'*char')'   ;
dum5          = fread(fid,1,'*float64')   ;
time          = fread(fid,1,'*float64')   ;
sfn           = fread(fid,1,'int32')      ;
dum6          = fread(fid,1,'*float64')   ;

aa=fread(fid,Inf,'*float64'); 

%Assign points taking into account gaps among processors
for i=1:ncores
    bb(1+sum(history2(1:i-1,2)):sum(history2(1:i,2)),1)=aa(1+sum(history2(1:i-1,2))+(i-1):sum(history2(1:i,2))+(i-1));
end

WD{field_number}=bb; 

fclose(fid);

%L2
%%%%%
fname=['ZSTAT/L',num2str(field_number,'%2.2d')];

[fid,message] = fopen(fname,'r','ieee-le');
hdr           = fread(fid,1,'int32')      ;
F             = fread(fid,hdr,'*char')'   ;
dum5          = fread(fid,1,'*float64')   ;
time          = fread(fid,1,'*float64')   ;
sfn           = fread(fid,1,'int32')      ;
dum6          = fread(fid,1,'*float64')   ;

aa=fread(fid,Inf,'*float64'); 

%Assign points taking into account gaps among processors
for i=1:ncores
    bb(1+sum(history2(1:i-1,2)):sum(history2(1:i,2)),1)=aa(1+sum(history2(1:i-1,2))+(i-1):sum(history2(1:i,2))+(i-1));
end

LD{field_number}=bb; 

fclose(fid);

%V
%%%%%
fname=['ZSTAT/V',num2str(field_number,'%2.2d')];

[fid,message] = fopen(fname,'r','ieee-le');
hdr           = fread(fid,1,'int32')      ;
F             = fread(fid,hdr,'*char')'   ;
dum5          = fread(fid,1,'*float64')   ;
time          = fread(fid,1,'*float64')   ;
sfn           = fread(fid,1,'int32')      ;
dum6          = fread(fid,1,'*float64')   ;

aa=fread(fid,Inf,'*float64'); 

%Assign points taking into account gaps among processors
for i=1:ncores
    bb(1+sum(history2(1:i-1,2)):sum(history2(1:i,2)),1)=aa(1+sum(history2(1:i-1,2))+(i-1):sum(history2(1:i,2))+(i-1));
end

VD{field_number}=bb; 

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure
% plot(x,U,'ob'); hold on;
% plot(x,UU,'or'); hold on;
% plot(x,UUU,'og'); hold on;
% 
% figure
% plot(x,V,'ob'); hold on;
% plot(x,VV,'or'); hold on;
% plot(x,VVV,'og'); hold on;
% 
% figure
% plot(x,W,'ob'); hold on;
% plot(x,WW,'or'); hold on;
% plot(x,WWW,'og'); hold on;
% 
% figure
% plot(phi,UUU,'ob'); hold on;
% 
% figure
% plot(phi,VVV,'or'); hold on;
% 
% figure
% plot(phi,WWW,'og'); hold on;

end

%% Define data
U = UD{1};
V = VD{1};
W = WD{1};
uu = LD{1};

vv = UD{2};
ww = VD{2};
uv = WD{2};
uw = LD{2};

vw = UD{3};
P = VD{3};
pp = WD{3};
ppp = LD{3};

pppp = UD{4};
uuu = VD{4};
vvv = WD{4};
www = LD{4};

uuv = UD{5};
uuw = VD{5};
uvv = WD{5};
vvw = LD{5};

uww = UD{6};
vww = VD{6};
uvw = WD{6};
Pxx = LD{6};

Pyy = UD{7};
Pzz = VD{7};
Pxy = WD{7};
Pxz = LD{7};

Pyz = UD{8};
Dxx = VD{8};
Dyy = WD{8};
Dzz = LD{8};

Dxy = UD{9};
Dxz = VD{9};
Dyz = WD{9};
Txx = LD{9};

Tyy = UD{10};
Tzz = VD{10};
Txy = WD{10};
Txz = LD{10};

Tyz = UD{11};
VDxx = VD{11};
VDyy = WD{11};
VDzz = LD{11};

VDxy = UD{12};
VDxz = VD{12};
VDyz = WD{12};
Pixx = LD{12};

Piyy = UD{13};
Pizz = VD{13};
Pixy = WD{13};
Pixz = LD{13};

Piyz = UD{14};
Cxx = VD{14};
Cyy = WD{14};
Czz = LD{14};

Cxy = UD{15};
Cxz = VD{15};
Cyz = WD{15};
Pk = LD{15};

Dk = UD{16};
Tk = VD{16};
VDk = WD{16};
Pik = LD{16};

Ck = UD{17};
Resk = VD{17};
PTxx = WD{17};
PTyy = LD{17};

PTzz = UD{18};
PTxy = VD{18};
PTxz = WD{18};
PTyz = LD{18};

PSxx = UD{19};
PSyy = VD{19};
PSzz = WD{19};
PSxy = LD{19};

PSxz = UD{20};
PSyz = VD{20};
dUdx = WD{20};
dUdy = LD{20};

dUdz = UD{21};
dVdx = VD{21};
dVdy = WD{21};
dVdz = LD{21};

dWdx = UD{22};
dWdy = VD{22};
dWdz = WD{22};

Re = 700;
nu = 1/Re;

save(['data_stats.mat'],'x','y','z',...
    'U','V','W','uu','vv','ww','uv','uw','vw','P','pp','ppp','pppp',...
    'uuu','vvv','www','uuv','uuw','uvv','vvw','uww','vww','uvw',...
    'Pxx','Pyy','Pzz','Pxy','Pxz','Pyz','Dxx','Dyy','Dzz','Dxy',...
    'Dxz','Dyz','Txx','Tyy','Tzz','Txy','Txz','Tyz','VDxx','VDyy',...
    'VDzz','VDxy','VDxz','VDyz','Pixx','Piyy','Pizz','Pixy','Pixz',...
    'Piyz','Cxx','Cyy','Czz','Cxy','Cxz','Cyz','Pk','Dk','Tk',...
    'VDk','Pik','Ck','Resk','PTxx','PTyy','PTzz','PTxy','PTxz',...
    'PTyz','PSxx','PSyy','PSzz','PSxy','PSxz','PSyz','dUdx','dUdy',...
    'dUdz','dVdx','dVdy','dVdz','dWdx','dWdy','dWdz','nu');










close all
clear all
clc
double precision;

%Initialize nponits
npoints=0;

%Domain
xmax = 9;
ymax = 3;

%% Generate interpolating mesh
x = linspace(0.1,xmax-0.1,40);
nmax = 42;
yn = 1-cos((2*(1:nmax)-1)/(2*nmax)*pi/2);
xx = zeros(length(x),length(yn));
yy = zeros(length(x),length(yn));

[y,alpha] = hill_values(x);

for i=1:length(x)
    xx(i,:) = x(i)+yn*sin(-alpha(i));
    yy(i,:) = y(i)+yn*cos(alpha(i));
end
x_pts = xx(:);
y_pts = yy(:);

npoints = length(x_pts)

%% Save data
%Save x data points in Fortran binary format
fid=fopen('./ZSTAT/x.fort','w','ieee-le.l64');

%First write 4 bytes integer
data=npoints;
eor=length(data)*4;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'int32');
count=fwrite(fid,eor,'int32');

%Then write npoints reals
data=x_pts;
eor=length(data)*8;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'float64');
count=fwrite(fid,eor,'int32');
fclose(fid);

%Save y data points in Fortran binary format
fid=fopen('./ZSTAT/y.fort','w','ieee-le.l64');

%First write 4 bytes integer
data=npoints;
eor=length(data)*4;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'int32');
count=fwrite(fid,eor,'int32');

%Then write npoints reals
data=y_pts;
eor=length(data)*8;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'float64');
count=fwrite(fid,eor,'int32');
fclose(fid);

figure(1)
plot(x_pts,y_pts,'.b'); hold on;

axis equal

save('int_mesh.mat');
close all; 
clear all; 
clc
double precision;
format long;

%Periodic hill geometry
xmin = 0;
xmax = 9;
ymin = 0;
ymax = 3.035;
zmin = 0;
zmax = 4.5;

%Wall normal positions
nmax = 30;
k = 1:nmax;
x = linspace(0.1,xmax-0.1,20);
yn = [0 1-cos((2*(1:nmax)-1)/(2*nmax)*pi/2) 1];
zn = linspace(zmin,zmax,10);

np = length(zn)*length(yn)*length(x)

xx = zeros(length(x),length(yn));
yy = zeros(length(x),length(yn));

[y,alpha] = hill_values(x);

for i=1:length(x)
    xx(i,:) = x(i)+yn*sin(-alpha(i));
    yy(i,:) = y(i)+yn*cos(alpha(i));
end
x_pts = xx(:);
y_pts = yy(:);


%Create coordinates
x=zeros(np,1);
y=zeros(np,1);
z=zeros(np,1);

for i = 1:length(zn)
    i
    f=length(x_pts);
    x(1+(i-1)*f:i*f) = x_pts;
    y(1+(i-1)*f:i*f) = y_pts;
    z(1+(i-1)*f:i*f) = zn(i);
end

figure
plot(x,y,'ob'); hold on;

figure
plot3(z,x,y,'ob'); hold on;
axis equal
xlabel('z')
ylabel('x')
zlabel('y')

%Save x data points in Fortran binary format
fid=fopen('ZSTAT/x.fort','w','ieee-le.l64');

%First write 4 bytes integer
data=np;
eor=length(data)*4;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'int32');
count=fwrite(fid,eor,'int32');

%Then write npoints reals
data=x;
eor=length(data)*8;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'float64');
count=fwrite(fid,eor,'int32');
fclose(fid);

%Save y data points in Fortran binary format
fid=fopen('ZSTAT/y.fort','w','ieee-le.l64');

%First write 4 bytes integer
data=np;
eor=length(data)*4;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'int32');
count=fwrite(fid,eor,'int32');

%Then write npoints reals
data=y;
eor=length(data)*8;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'float64');
count=fwrite(fid,eor,'int32');
fclose(fid);

%Save z data points in Fortran binary format
fid=fopen('ZSTAT/z.fort','w','ieee-le.l64');

%First write 4 bytes integer
data=np;
eor=length(data)*4;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'int32');
count=fwrite(fid,eor,'int32');

%Then write npoints reals
data=z;
eor=length(data)*8;
count=fwrite(fid,eor,'int32');
count=fwrite(fid,data,'float64');
count=fwrite(fid,eor,'int32');
fclose(fid);

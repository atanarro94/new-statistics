%close all; 
clear all; 
clc
double precision;
format long;

ifplt = 0

% Fields in the binary record from stat files

%Statistic fields
% u,v,w,p are instantaneous quantities 
% averaged in time and homogeneous direction z 

% 1.  <u>           % F1
% 2.  <v>           % F2
% 3.  <w>           % F3
% 4.  <p>           % F4

% 5.  <uu>          % F5
% 6.  <vv>          % F6
% 7.  <ww>          % F7
% 8.  <pp>          % F8

% 9.  <uv>          % F9
% 10. <vw>          % F10
% 11. <uw>          % F11

% 12. <pu>          % F12
% 13. <pv>          % F13
% 14. <pw>          % F14

% 15. <pdudx>       % F15
% 16. <pdudy>       % F16
% 17. <pdudz>       % F17

% 18. <pdvdx>       % F18
% 19. <pdvdy>       % F19
% 20. <pdvdz>       % F20

% 21. <pdwdx>       % F21
% 22. <pdwdy>       % F22
% 23. <pdwdz>       % F23

% 24. <uuu>         % F24
% 25. <vvv>         % F25
% 26. <www>         % F26
% 27. <ppp>         % F27          

% 28. <uuv>         % F28          
% 29. <uuw>         % F29          
% 30. <vvu>         % F30          
% 31. <vvw>  	     % F31             
% 32. <wwu>         % F32          
% 33. <wwv>         % F33          
% 34. <uvw>         % F34          

% 35. <uuuu>        % F35          
% 36. <vvvv>        % F36          
% 37. <wwww>        % F37          
% 38. <pppp>        % F38          

% 39. e11: <((du/dx)^2+(du/dy)^2+(du/dz)^2)>                % F39
% 40. e22: <((dv/dx)^2+(dv/dy)^2+(dv/dz)^2)>                % F40 
% 41. e33: <((dw/dx)^2+(dw/dy)^2+(dw/dz)^2)>                % F41 
% 42. e12: <(du/dx*dv/dx+du/dy*dv/dy+du/dz*dv/dz)>          % F42  
% 43. e13: <(du/dx*dw/dx+du/dy*dw/dy+du/dz*dw/dz)>          % F43 
% 44. e23: <(dv/dx*dw/dx+dv/dy*dw/dy+dv/dz*dw/dz)>          % F44 

%Derivative fields
% 1. dU/dx           % D1
% 2. dU/dy           % D2
% 3. dV/dx           % D3
% 4. dV/dy           % D4

% 5. dW/dx           % D5
% 6. dW/dy           % D6
% 7. dP/dx           % D7
% 8. dP/dy           % D8

% 9.  d<uu>/dx       % D9
% 10. d<uu>/dy       % D10
% 11. d<vv>/dx       % D11
% 12. d<vv>/dy       % D12

% 13. d<ww>/dx       % D13
% 14. d<ww>/dy       % D14
% 15. d<pp>/dx       % D15
% 16. d<pp>/dy       % D16

% 17. d<uv>/dx       % D17
% 18. d<uv>/dy       % D18
% 19. d<vw>/dx       % D19
% 20. d<vw>/dy       % D20

% 21. d<uw>/dx       % D21
% 22. d<uw>/dy       % D22
% 23. d<uuu>/dx      % D23
% 24. d<uuu>/dy      % D24

% 25. d<vvv>/dx      % D25
% 26. d<vvv>/dy      % D26
% 27. d<www>/dx      % D27
% 28. d<www>/dy      % D28

% 29. d<ppp>/dx      % D29
% 30. d<ppp>/dy      % D30
% 31. d<uuv>/dx      % D31
% 32. d<uuv>/dy      % D32

% 33. d<uuw>/dx      % D33
% 34. d<uuw>/dy      % D34
% 35. d<vvu>/dx      % D35
% 36. d<vvu>/dy      % D36

% 37. d<vvw>/dx      % D37
% 38. d<vvw>/dy      % D38
% 39. d<wwu>/dx      % D39
% 40. d<wwu>/dy      % D40

% 41. d<wwv>/dx      % D41
% 42. d<wwv>/dy      % D42
% 43. d<uvw>/dx      % D43
% 44. d<uvw>/dy      % D44

% 45. d2U/dx2        % D45
% 46. d2U/dy2        % D46
% 47. d2V/dx2        % D47
% 48. d2V/dy2        % D48

% 49. d2W/dx2        % D49
% 50. d2W/dy2        % D50
% 51. d2<uu>/dx2     % D51
% 52. d2<uu>/dy2     % D52

% 53. d2<vv>/dx2     % D53
% 54. d2<vv>/dy2     % D54
% 55. d2<ww>/dx2     % D55
% 56. d2<ww>/dy2     % D56

% 57. d2<uv>/dx2     % D57
% 58. d2<uv>/dy2     % D58
% 59. d2<uw>/dx2     % D59
% 60. d2<uw>/dy2     % D60

% 61. d2<vw>/dx2     % D61
% 62. d2<vw>/dy2     % D62

% 63. d<pu>/dx       % D63
% 64. d<pu>/dy       % D64
% 65. d<pv>/dx       % D65
% 66. d<pv>/dy       % D66

% 67. d<pw>/dx       % D67
% 68. d<pw>/dy       % D68

% Read interpolated field: int_fld file
folder        = '/scratch/atdr/simulations/periodic_hill_old_stats/post2d/';
fsave         = 'data_old_stat.mat';
fname         = [folder,'ZSTAT/int_fld'];
[fid,message] = fopen(fname,'r','ieee-le');
hdr           = fread(fid,1,'int32')      ;
F             = fread(fid,hdr,'*char')'   ;
dum5          = fread(fid,1,'*float64')   ;
Rer           = fread(fid,1,'*float64')   ;
Domain        = fread(fid,3,'*float64')   ;
nel           = fread(fid,3,'int32')      ;
Poly          = fread(fid,3,'int32')      ;
nstat         = fread(fid,1,'int32')      ;
nderiv        = fread(fid,1,'int32')      ;
times         = fread(fid,1,'*float64')   ;
timee         = fread(fid,1,'*float64')   ;
atime         = fread(fid,1,'*float64')   ;
DT            = fread(fid,1,'*float64')   ;
nrec          = fread(fid,1,'int32')      ;
Tint          = fread(fid,1,'*float64')   ;
npoints       = fread(fid,1,'int32')      ;
dum6          = fread(fid,1,'*float64')   ;

% Read interpolating mesh: x.fort and y.fort
fname_x           = [folder,'ZSTAT/x.fort'];
[fid_x,message_x] = fopen(fname_x,'r','ieee-le');
hdr_x             = fread(fid_x,1,'int32');
dum1              = fread(fid_x,3,'int32');
% dum2              = fread(fid_x,1,'*float64');
x_pts             = fread(fid_x,npoints,'*float64');
fclose(fid_x);

fname_y          = [folder,'ZSTAT/y.fort'];
[fid_y,message_y] = fopen(fname_y,'r','ieee-le');
hdr_y             = fread(fid_y,1,'int32');
dum3              = fread(fid_y,3,'int32');
% dum4              = fread(fid_y,1,'*float64');
y_pts             = fread(fid_y,npoints,'*float64');
fclose(fid_y);

% Define data matrices

nx=40;
ny=42;

X = reshape(x_pts,[nx ny]);
Y = reshape(y_pts,[nx ny]);
Mat=zeros(nx,ny);
D=F;

% This is how data is stored in matrix form
%
%          x=-1                                      x=+5.2 
% y = -1 [ |*|     |*|     |*|     |*|     |*|     |*|   ] 
%        [ |*|      *       *       *       *       *    ] 
%        [ |*|      *       *       *       *       *    ]
%        [ |*|      *       *       *       *       *    ]
% y = +1 [ |*|      *       *       *       *       *    ] 


% Arrange stat fields in matrix form
for ii=1:nstat+nderiv
    
    if ii == 1
        fseek(fid,0,'cof');
    else
        fseek(fid,8,'cof');
    end

    % Store current field as an array
    MatR = fread(fid,npoints,'*float64'); 

    % Arrange current field in matrix form
    for i=1:ny
        Mat(:,i)=MatR((i-1)*nx+1:i*nx);
    end
    
    if ii<=nstat
        % Generata variable name for field X as FX
        v=genvarname('F', who);
    else
        % Generata variable name for field X as DX
        v=genvarname('D', who);
    end
    
    % Store matrix in variable FX or Dx
    evalc([v '=INOUT(Mat)']);   
    
end

fclose(fid);

% Simulation parameters
% Bulk Reynolds number Reb=Ub*c/nu
Reb=Rer; 

% Domain dimensions
Lx=Domain(1);
Ly=Domain(2);

% Fluid density
rho=1; 

% Kinematic viscosity. Both Ub and c are unit normalizing parameters
nu=1/Reb;

% Molecular visosity
mu=rho*nu;

%Mean velocities. Tensors of Rank 1.
R1_tensor(1:3,1,1:nx,1:ny)=0;

for j=1:ny
    for i=1:nx
        R1_tensor(:,1,i,j)=[F1(i,j); F2(i,j); F3(i,j)];
    end
end

U=squeeze(squeeze(R1_tensor(1,1,:,:)));
V=squeeze(squeeze(R1_tensor(2,1,:,:)));
W=squeeze(squeeze(R1_tensor(3,1,:,:)));

%Reynolds stress tensor. Tensor of Rank 2.
R2_tensor(1:3,1:3,1:nx,1:ny)=0;
for j=1:ny
    for i=1:nx
        R2_tensor(:,:,i,j)=[ (F5(i,j)-F1(i,j).*F1(i,j)) (F9(i,j)-F1(i,j).*F2(i,j)) (F11(i,j)-F1(i,j).*F3(i,j));
                             (F9(i,j)-F1(i,j).*F2(i,j)) (F6(i,j)-F2(i,j).*F2(i,j)) (F10(i,j)-F2(i,j).*F3(i,j)); 
                            (F11(i,j)-F1(i,j).*F3(i,j)) (F10(i,j)-F2(i,j).*F3(i,j)) (F7(i,j)-F3(i,j).*F3(i,j))];
    end    
end

uu=squeeze(squeeze(R2_tensor(1,1,:,:)));
vv=squeeze(squeeze(R2_tensor(2,2,:,:)));
ww=squeeze(squeeze(R2_tensor(3,3,:,:)));
uv=squeeze(squeeze(R2_tensor(1,2,:,:)));
uw=squeeze(squeeze(R2_tensor(1,3,:,:)));
vw=squeeze(squeeze(R2_tensor(2,3,:,:)));

%Mean, RMS, skewness and flatness of pressure
P=F4;
pp=F8-P.*P;
ppp=F27-3*P.*pp-P.*P.*P;
pppp=F38-4*P.*ppp-6*P.*P.*pp-P.*P.*P.*P;

%Normalize pressure
prms=sqrt(pp);
pskew=ppp./(pp).^(3/2);
pflat=pppp./(pp).^(2);

%Skewness tensor. Tensor of Rank 3.
R3_tensor_tot(1:3,1:3,1:3,1:nx,1:ny)=0;
R3_tensor(1:3,1:3,1:3,1:nx,1:ny)=0;

%Form of the tensor.
% [ uuu uuv uuw ] [ uuv uvv uvw ] [ uuw uvw uww ] 
% [ uuv uvv uvw ] [ uvv vvv vvw ] [ uvw vvw vww ]
% [ uuw uvw uww ] [ uvw vvw vww ] [ uww vww www ]

for j=1:ny
    for i=1:nx
        R3_tensor_tot(:,:,1,i,j) = [F24(i,j) F28(i,j) F29(i,j); 
                                    F28(i,j) F30(i,j) F34(i,j); 
                                    F29(i,j) F34(i,j) F32(i,j)];
             
        R3_tensor_tot(:,:,2,i,j) = [F28(i,j) F30(i,j) F34(i,j); 
                                    F30(i,j) F25(i,j) F31(i,j); 
                                    F34(i,j) F31(i,j) F33(i,j)];
             
        R3_tensor_tot(:,:,3,i,j) = [F29(i,j) F34(i,j) F32(i,j); 
                                    F34(i,j) F31(i,j) F33(i,j); 
                                    F32(i,j) F33(i,j) F26(i,j)]; 
    end    
end

for j=1:ny
    for i=1:nx
        R3_tensor(:,:,1,i,j)=[(R3_tensor_tot(1,1,1,i,j)-3*U(i,j).*uu(i,j)-U(i,j).*U(i,j).*U(i,j)) ...
                              (R3_tensor_tot(1,2,1,i,j)-2*U(i,j).*uv(i,j)-V(i,j).*uu(i,j)-U(i,j).*U(i,j).*V(i,j)) ... 
                              (R3_tensor_tot(1,3,1,i,j)-2*U(i,j).*uw(i,j)-W(i,j).*uu(i,j)-U(i,j).*U(i,j).*W(i,j));
                              (R3_tensor_tot(1,2,1,i,j)-2*U(i,j).*uv(i,j)-V(i,j).*uu(i,j)-U(i,j).*U(i,j).*V(i,j)) ... 
                              (R3_tensor_tot(2,2,1,i,j)-2*V(i,j).*uv(i,j)-U(i,j).*vv(i,j)-V(i,j).*V(i,j).*U(i,j)) ...
                              (R3_tensor_tot(2,3,1,i,j)-U(i,j).*vw(i,j)-V(i,j).*uw(i,j)-W(i,j).*uv(i,j)-U(i,j).*V(i,j).*W(i,j));
                              (R3_tensor_tot(1,3,1,i,j)-2*U(i,j).*uw(i,j)-W(i,j).*uu(i,j)-U(i,j).*U(i,j).*W(i,j)) ...
                              (R3_tensor_tot(2,3,1,i,j)-U(i,j).*vw(i,j)-V(i,j).*uw(i,j)-W(i,j).*uv(i,j)-U(i,j).*V(i,j).*W(i,j)) ...
                              (R3_tensor_tot(3,3,1,i,j)-2*W(i,j).*uw(i,j)-U(i,j).*ww(i,j)-W(i,j).*W(i,j).*U(i,j))];
                           
        R3_tensor(:,:,2,i,j)=[(R3_tensor_tot(1,1,2,i,j)-2*U(i,j).*uv(i,j)-V(i,j).*uu(i,j)-U(i,j).*U(i,j).*V(i,j)) ... 
                              (R3_tensor_tot(1,2,2,i,j)-2*V(i,j).*uv(i,j)-U(i,j).*vv(i,j)-V(i,j).*V(i,j).*U(i,j)) ...
                              (R3_tensor_tot(1,3,2,i,j)-U(i,j).*vw(i,j)-V(i,j).*uw(i,j)-W(i,j).*uv(i,j)-U(i,j).*V(i,j).*W(i,j));
                              (R3_tensor_tot(1,2,2,i,j)-2*V(i,j).*uv(i,j)-U(i,j).*vv(i,j)-V(i,j).*V(i,j).*U(i,j)) ...
                              (R3_tensor_tot(2,2,2,i,j)-3*V(i,j).*vv(i,j)-V(i,j).*V(i,j).*V(i,j)) ...
                              (R3_tensor_tot(2,3,2,i,j)-2*V(i,j).*vw(i,j)-W(i,j).*vv(i,j)-V(i,j).*V(i,j).*W(i,j));
                              (R3_tensor_tot(1,3,2,i,j)-U(i,j).*vw(i,j)-V(i,j).*uw(i,j)-W(i,j).*uv(i,j)-U(i,j).*V(i,j).*W(i,j)) ...
                              (R3_tensor_tot(2,3,2,i,j)-2*V(i,j).*vw(i,j)-W(i,j).*vv(i,j)-V(i,j).*V(i,j).*W(i,j)) ...
                              (R3_tensor_tot(3,3,2,i,j)-2*W(i,j).*vw(i,j)-V(i,j).*ww(i,j)-W(i,j).*W(i,j).*V(i,j))];
        
                           
        R3_tensor(:,:,3,i,j)=[(R3_tensor_tot(1,1,3,i,j)-2*U(i,j).*uw(i,j)-W(i,j).*uu(i,j)-U(i,j).*U(i,j).*W(i,j)) ...
                              (R3_tensor_tot(1,2,3,i,j)-U(i,j).*vw(i,j)-V(i,j).*uw(i,j)-W(i,j).*uv(i,j)-U(i,j).*V(i,j).*W(i,j)) ...
                              (R3_tensor_tot(1,3,3,i,j)-2*W(i,j).*uw(i,j)-U(i,j).*ww(i,j)-W(i,j).*W(i,j).*U(i,j));
                              (R3_tensor_tot(1,2,3,i,j)-U(i,j).*vw(i,j)-V(i,j).*uw(i,j)-W(i,j).*uv(i,j)-U(i,j).*V(i,j).*W(i,j)) ...
                              (R3_tensor_tot(2,2,3,i,j)-2*V(i,j).*vw(i,j)-W(i,j).*vv(i,j)-V(i,j).*V(i,j).*W(i,j)) ...
                              (R3_tensor_tot(2,3,3,i,j)-2*W(i,j).*vw(i,j)-V(i,j).*ww(i,j)-W(i,j).*W(i,j).*V(i,j));
                              (R3_tensor_tot(1,3,3,i,j)-2*W(i,j).*uw(i,j)-U(i,j).*ww(i,j)-W(i,j).*W(i,j).*U(i,j)) ...
                              (R3_tensor_tot(2,3,3,i,j)-2*W(i,j).*vw(i,j)-V(i,j).*ww(i,j)-W(i,j).*W(i,j).*V(i,j)) ...
                              (R3_tensor_tot(3,3,3,i,j)-3*W(i,j).*ww(i,j)-W(i,j).*W(i,j).*W(i,j))];
    end
end

uuu=squeeze(squeeze(R3_tensor(1,1,1,:,:)));
vvv=squeeze(squeeze(R3_tensor(2,2,2,:,:)));
www=squeeze(squeeze(R3_tensor(3,3,3,:,:)));
uuv=squeeze(squeeze(R3_tensor(1,2,1,:,:)));
uuw=squeeze(squeeze(R3_tensor(1,3,1,:,:)));
uvv=squeeze(squeeze(R3_tensor(2,2,1,:,:)));
vvw=squeeze(squeeze(R3_tensor(2,3,2,:,:)));
uww=squeeze(squeeze(R3_tensor(3,3,1,:,:)));
vww=squeeze(squeeze(R3_tensor(3,3,2,:,:)));
uvw=squeeze(squeeze(R3_tensor(2,3,1,:,:)));

%Normalized skewness
uskew=uuu./real(uu.^(3/2));
vskew=vvv./real(vv.^(3/2));
wskew=www./real(ww.^(3/2));

%Velocity gradient tensor. Tensor of Rank 2.
dUidxj(1:3,1:3,1:nx,1:ny)=0;

dUdx=D1;
dVdx=D3;
dWdx=D5;

dUdy=D2;
dVdy=D4;
dWdy=D6;

dUdz(1:nx,1:ny)=0;
dVdz(1:nx,1:ny)=0;
dWdz(1:nx,1:ny)=0;

for j=1:ny
    for i=1:nx
        dUidxj(:,:,i,j)=[dUdx(i,j) dUdy(i,j) dUdz(i,j);
                         dVdx(i,j) dVdy(i,j) dVdz(i,j);
                         dWdx(i,j) dWdy(i,j) dWdz(i,j)];
    end
end


%Production tensor. Tensor of Rank 2.
P_tensor(1:3,1:3,1:nx,1:ny)=0;

%Definition of the production tensor assuming fully-developed flow, i.e.,
%d()/dz=0, as a function of x and y.
%Pij=-(<uiuk>*dUj/dxk+<ujuk>*dUi/dxk)
%P11=-2*(<uu>*dU/dx+<uv>*dU/dy) 
%P12=-(<uu>*dV/dx +<uv>*dV/dy+<uv>*dU/dx+<vv>*dU/dy)
%P13=-(<uu>*dW/dx +<uv>*dW/dy+<uw>*dU/dx+<vw>*dU/dy)
%P22=-2*(<uv>*dV/dx+<vv>*dV/dy) 
%P23=-(<uv>*dW/dx+<vv>*dW/dy+<uw>*dV/dx+<vw>*dV/dy)
%P33=-2*(<uw>*dW/dx+<vw>*dW/dy)

for j=1:ny
    for i=1:nx
        P_tensor(:,:,i,j)=[-2*(uu(i,j).*dUdx(i,j)+uv(i,j).*dUdy(i,j)) ...
                           -(uu(i,j).*dVdx(i,j)+uv(i,j).*dVdy(i,j)+uv(i,j).*dUdx(i,j)+vv(i,j).*dUdy(i,j)) ...
                           -(uu(i,j).*dWdx(i,j)+uv(i,j).*dWdy(i,j)+uw(i,j).*dUdx(i,j)+vw(i,j).*dUdy(i,j));
                           -(uu(i,j).*dVdx(i,j)+uv(i,j).*dVdy(i,j)+uv(i,j).*dUdx(i,j)+vv(i,j).*dUdy(i,j)) ...
                           -2*(uv(i,j).*dVdx(i,j)+vv(i,j).*dVdy(i,j)) ...
                           -(uv(i,j).*dWdx(i,j)+vv(i,j).*dWdy(i,j)+uw(i,j).*dVdx(i,j)+vw(i,j).*dVdy(i,j));
                           -(uu(i,j).*dWdx(i,j)+uv(i,j).*dWdy(i,j)+uw(i,j).*dUdx(i,j)+vw(i,j).*dUdy(i,j)) ...
                           -(uv(i,j).*dWdx(i,j)+vv(i,j).*dWdy(i,j)+uw(i,j).*dVdx(i,j)+vw(i,j).*dVdy(i,j)) ...
                           -2*(uw(i,j).*dWdx(i,j)+vw(i,j).*dWdy(i,j))];
    end
end

Pxx=squeeze(squeeze(P_tensor(1,1,:,:)));
Pyy=squeeze(squeeze(P_tensor(2,2,:,:)));
Pzz=squeeze(squeeze(P_tensor(3,3,:,:)));
Pxy=squeeze(squeeze(P_tensor(1,2,:,:)));
Pyz=squeeze(squeeze(P_tensor(2,3,:,:)));
Pxz=squeeze(squeeze(P_tensor(1,3,:,:)));

%Dissipation tensor. Tensor of Rank 2.
D_tensor(1:3,1:3,1:nx,1:ny)=0;

%Definition of the dissipation tensor assuming fully-developed flow, i.e.,
%d()/dz=0, as a function of x and y.
%Dij=-2*nu*<dui/dxk*duj/dxk>

%e11_tot=<((du/dx)^2+(du/dy)^2+(du/dz)^2)> % F39 
%e22_tot=<((dv/dx)^2+(dv/dy)^2+(dv/dz)^2)> % F40 
%e33_tot=<((dw/dx)^2+(dw/dy)^2+(dw/dz)^2)> % F41 
%e12_tot=<(du/dx*dv/dx+du/dy*dv/dy+du/dz*dv/dz)>  % F42  
%e13_tot=<(du/dx*dw/dx+du/dy*dw/dy+du/dz*dw/dz)>  % F43 
%e23_tot=<(dv/dx*dw/dx+dv/dy*dw/dy+dv/dz*dw/dz)>  % F44 

%e11=e11_tot-(dU/dx)^2-(dU/dy)^2 
%e22=e22_tot-(dV/dx)^2-(dV/dy)^2 
%e33=e33_tot-(dW/dx)^2-(dW/dy)^2 
%e12=e12_tot-(dU/dx*dV/x)-(dU/dy*dV/dy)
%e13=e13_tot-(dU/dx*dW/x)-(dU/dy*dW/dy)
%e23=e23_tot-(dV/dx*dW/x)-(dV/dy*dW/dy)

e11(1:nx,1:ny)=0;
e22(1:nx,1:ny)=0;
e33(1:nx,1:ny)=0;
e12(1:nx,1:ny)=0;
e13(1:nx,1:ny)=0;
e23(1:nx,1:ny)=0;

for j=1:ny
    e11(:,j)=F39(:,j)- ...
        (squeeze(dUidxj(1,1,:,j))).^2-(squeeze(dUidxj(1,2,:,j))).^2; 
    e22(:,j)=F40(:,j)- ...
        (squeeze(dUidxj(2,1,:,j))).^2-(squeeze(dUidxj(2,2,:,j))).^2; 
    e33(:,j)=F41(:,j)- ...
        (squeeze(dUidxj(3,1,:,j))).^2-(squeeze(dUidxj(3,2,:,j))).^2; 
    e12(:,j)=F42(:,j)- ...
        squeeze(dUidxj(1,1,:,j)).*squeeze(dUidxj(2,1,:,j))- ...
        squeeze(dUidxj(1,2,:,j)).*squeeze(dUidxj(2,2,:,j));
    e13(:,j)=F43(:,j)- ...
        squeeze(dUidxj(1,1,:,j)).*squeeze(dUidxj(3,1,:,j))- ...
        squeeze(dUidxj(1,2,:,j)).*squeeze(dUidxj(3,2,:,j));
    e23(:,j)=F44(:,j)- ...
        squeeze(dUidxj(2,1,:,j)).*squeeze(dUidxj(3,1,:,j))- ...
        squeeze(dUidxj(2,2,:,j)).*squeeze(dUidxj(3,2,:,j));      
end

for j=1:ny
    for i=1:nx
        D_tensor(:,:,i,j)=[e11(i,j) e12(i,j) e13(i,j);
                           e12(i,j) e22(i,j) e23(i,j);
                           e13(i,j) e23(i,j) e33(i,j)];
    end
end

Dxx=-2*nu*squeeze(squeeze(D_tensor(1,1,:,:)));
Dyy=-2*nu*squeeze(squeeze(D_tensor(2,2,:,:)));
Dzz=-2*nu*squeeze(squeeze(D_tensor(3,3,:,:)));
Dxy=-2*nu*squeeze(squeeze(D_tensor(1,2,:,:)));
Dyz=-2*nu*squeeze(squeeze(D_tensor(2,3,:,:)));
Dxz=-2*nu*squeeze(squeeze(D_tensor(1,3,:,:)));

%Mean convection tensor. Tensor of Rank 2.
C_tensor(1:3,1:3,1:nx,1:ny)=0;

%Definition of the mean convection tensor assuming 
%fully-developed flow, i.e., d()/dz=0, as a function of x and y.
%Cij=Uk*d<uiuj>/dxk
%Note that under this definition: Production + Dissipation - Convection ...
%C11=U*d(uu)/dx+V*d(uu)/dy
%C22=U*d(vv)/dx+V*d(vv)/dy
%C33=U*d(ww)/dx+V*d(ww)/dy
%C12=U*d(uv)/dx+V*d(uv)/dy
%C13=U*d(uw)/dx+V*d(uw)/dy
%C23=U*d(vw)/dx+V*d(vw)/dy

duudx=D9-2*U.*dUdx;
dvvdx=D11-2*V.*dVdx;
dwwdx=D13-2*W.*dWdx;
duvdx=D17-U.*dVdx-V.*dUdx;
duwdx=D21-U.*dWdx-W.*dUdx;
dvwdx=D19-V.*dWdx-W.*dVdx;

duudy=D10-2*U.*dUdy;
dvvdy=D12-2*V.*dVdy;
dwwdy=D14-2*W.*dWdy;
duvdy=D18-U.*dVdy-V.*dUdy;
duwdy=D22-U.*dWdy-W.*dUdy;
dvwdy=D20-V.*dWdy-W.*dVdy;

duudz(1:nx,1:ny)=0;
dvvdz(1:nx,1:ny)=0;
dwwdz(1:nx,1:ny)=0;
duvdz(1:nx,1:ny)=0;
duwdz(1:nx,1:ny)=0;
dvwdz(1:nx,1:ny)=0;

for j=1:ny
    for i=1:nx
        C_tensor(:,:,i,j)=[U(i,j).*duudx(i,j)+V(i,j).*duudy(i,j) ...
                           U(i,j).*duvdx(i,j)+V(i,j).*duvdy(i,j) ...
                           U(i,j).*duwdx(i,j)+V(i,j).*duwdy(i,j);
                           U(i,j).*duvdx(i,j)+V(i,j).*duvdy(i,j) ...
                           U(i,j).*dvvdx(i,j)+V(i,j).*dvvdy(i,j) ...
                           U(i,j).*dvwdx(i,j)+V(i,j).*dvwdy(i,j);
                           U(i,j).*duwdx(i,j)+V(i,j).*duwdy(i,j) ...
                           U(i,j).*dvwdx(i,j)+V(i,j).*dvwdy(i,j) ...
                           U(i,j).*dwwdx(i,j)+V(i,j).*dwwdy(i,j)];
    end
end

Cxx=squeeze(squeeze(C_tensor(1,1,:,:)));
Cyy=squeeze(squeeze(C_tensor(2,2,:,:)));
Czz=squeeze(squeeze(C_tensor(3,3,:,:)));
Cxy=squeeze(squeeze(C_tensor(1,2,:,:)));
Cxz=squeeze(squeeze(C_tensor(1,3,:,:)));
Cyz=squeeze(squeeze(C_tensor(2,3,:,:)));

%Turbulent transport tensor. Tensor of Rank 2.
T_tensor(1:3,1:3,1:nx,1:ny)=0;

%Definition of the turbulent transport tensor assuming fully-developed 
%flow, i.e., d()/dz=0, as a function of x and y.
%Tij=-d<uiujuk>/dxk
%T11=-(d<uuuk>)/dxk, K=1,3=-[d<uuu>/dx+d<uuv>/dy+d<uuw>/dz]
%T22=-(d<vvuk>)/dxk, K=1,3=-[d<vvu>/dx+d<vvv>/dy+d<vvw>/dz]
%T33=-(d<wwuk>)/dxk, K=1,3=-[d<wwu>/dx+d<wwv>/dy+d<www>/dz]
%T12=-(d<uvuk>)/dxk, K=1,3=-[d<uvu>/dx+d<uvv>/dy+d<uvw>/dz]
%T13=-(d<uwuk>)/dxk, K=1,3=-[d<uwu>/dx+d<uwv>/dy+d<uww>/dz]
%T23=-(d<vwuk>)/dxk, K=1,3=-[d<vwu>/dx+d<vwv>/dy+d<vww>/dz]

duuudx=D23-3*U.*U.*dUdx-3*(U.*duudx+uu.*dUdx);
dvvudx=D35-2*(V.*duvdx+uv.*dVdx)-(U.*dvvdx+vv.*dUdx) ...
    -(V.*V.*dUdx+2*U.*V.*dVdx);
dwwudx=D39-2*(W.*duwdx+uw.*dWdx)-(U.*dwwdx+ww.*dUdx) ... 
    -(W.*W.*dUdx+2*U.*W.*dWdx);
duvudx=D31-2*(U.*duvdx+uv.*dUdx)-(V.*duudx+uu.*dVdx) ...
    -(U.*U.*dVdx+2*U.*V.*dUdx);
duwudx=D33-2*(U.*duwdx+uw.*dUdx)-(W.*duudx+uu.*dWdx) ...
    -(U.*U.*dWdx+2*U.*W.*dUdx);
dvwudx=D43-(U.*dvwdx+vw.*dUdx)-(V.*duwdx+uw.*dVdx)-(W.*duvdx+uv.*dWdx) ...
    -(U.*V.*dWdx+U.*W.*dVdx+V.*W.*dUdx);
      
duuvdy=D32-2*(U.*duvdy+uv.*dUdy)-(V.*duudy+uu.*dVdy) ...
    -(U.*U.*dVdy+2*U.*V.*dUdy);
dvvvdy=D26-3*(V.*dvvdy+vv.*dVdy)-3*V.*V.*dVdy;
dwwvdy=D42-2*(W.*dvwdy+vw.*dWdy)-(V.*dwwdy+ww.*dVdy) ...
    -(W.*W.*dVdy+2*V.*W.*dWdy);
duvvdy=D36-2*(V.*duvdy+uv.*dVdy)-(U.*dvvdy+vv.*dUdy) ...
    -(V.*V.*dUdy+2*U.*V.*dVdy);
duwvdy=D44-(U.*dvwdy+vw.*dUdy)-(V.*duwdy+uw.*dVdy)-(W.*duvdy+uv.*dWdy) ...
    -(U.*V.*dWdy+U.*W.*dVdy+V.*W.*dUdy);
dvwvdy=D38-2*(V.*dvwdy+vw.*dVdy)-(W.*dvvdy+vv.*dWdy) ...
    -(V.*V.*dWdy+2*V.*W.*dVdy);

duuwdz(1:nx,1:ny)=0;
dvvwdz(1:nx,1:ny)=0;
dwwwdz(1:nx,1:ny)=0;
duvwdz(1:nx,1:ny)=0;
duwwdz(1:nx,1:ny)=0;
dvwwdz(1:nx,1:ny)=0;

for j=1:ny
    for i=1:nx
        T_tensor(:,:,i,j)=[-(duuudx(i,j)+duuvdy(i,j)+duuwdz(i,j)) ...
                           -(duvudx(i,j)+duvvdy(i,j)+duvwdz(i,j)) ...
                           -(duwudx(i,j)+duwvdy(i,j)+duwwdz(i,j));
                           -(duvudx(i,j)+duvvdy(i,j)+duvwdz(i,j)) ...
                           -(dvvudx(i,j)+dvvvdy(i,j)+dvvwdz(i,j)) ...
                           -(dvwudx(i,j)+dvwvdy(i,j)+dvwwdz(i,j));
                           -(duwudx(i,j)+duwvdy(i,j)+duwwdz(i,j)) ...
                           -(dvwudx(i,j)+dvwvdy(i,j)+dvwwdz(i,j)) ...
                           -(dwwudx(i,j)+dwwvdy(i,j)+dwwwdz(i,j))];
    end
end

Txx=squeeze(squeeze(T_tensor(1,1,:,:)));
Tyy=squeeze(squeeze(T_tensor(2,2,:,:)));
Tzz=squeeze(squeeze(T_tensor(3,3,:,:)));
Txy=squeeze(squeeze(T_tensor(1,2,:,:)));
Txz=squeeze(squeeze(T_tensor(1,3,:,:)));
Tyz=squeeze(squeeze(T_tensor(2,3,:,:)));

%Viscous diffusion tensor. Tensor of Rank 2.
VD_tensor(1:3,1:3,1:nx,1:ny)=0;

%Definition of the viscous diffusion tensor assuming fully-developed 
%flow, i.e., d()/dz=0, as a function of x and y.
%VDij=nu*d2(uiuj)/dxk2
%VD11=nu*(d2(uu)/dx2+d2(uu)/dy2+d2(uu)/dz2)
%VD22=nu*(d2(vv)/dx2+d2(vv)/dy2+d2(vv)/dz2)
%VD33=nu*(d2(ww)/dx2+d2(ww)/dy2+d2(ww)/dz2)
%VD12=nu*(d2(uv)/dx2+d2(uv)/dy2+d2(uv)/dz2)
%VD13=nu*(d2(uw)/dx2+d2(uw)/dy2+d2(uw)/dz2)
%VD23=nu*(d2(vw)/dx2+d2(vw)/dy2+d2(vw)/dz2)

d2uudx2=D51-2*(U.*D45+dUdx.*dUdx);
d2vvdx2=D53-2*(V.*D47+dVdx.*dVdx);
d2wwdx2=D55-2*(W.*D49+dWdx.*dWdx);
d2uvdx2=D57-(V.*D45+U.*D47+2*dUdx.*dVdx);
d2uwdx2=D59-(U.*D49+W.*D45+2*dUdx.*dWdx);
d2vwdx2=D61-(V.*D49+W.*D47+2*dVdx.*dWdx);

d2uudy2=D52-2*(U.*D46+dUdy.*dUdy);
d2vvdy2=D54-2*(V.*D48+dVdy.*dVdy);
d2wwdy2=D56-2*(W.*D50+dWdy.*dWdy);
d2uvdy2=D58-(V.*D46+U.*D48+2*dUdy.*dVdy);
d2uwdy2=D60-(U.*D50+W.*D46+2*dUdy.*dWdy);
d2vwdy2=D62-(V.*D50+W.*D48+2*dVdy.*dWdy);

d2uudz2(1:nx,1:ny)=0;
d2vvdz2(1:nx,1:ny)=0;
d2wwdz2(1:nx,1:ny)=0;
d2uvdz2(1:nx,1:ny)=0;
d2uwdz2(1:nx,1:ny)=0;
d2vwdz2(1:nx,1:ny)=0;

for j=1:ny
    for i=1:nx
        VD_tensor(:,:,i,j)=[nu*(d2uudx2(i,j)+d2uudy2(i,j)+d2uudz2(i,j)) ...
                            nu*(d2uvdx2(i,j)+d2uvdy2(i,j)+d2uvdz2(i,j)) ...
                            nu*(d2uwdx2(i,j)+d2uwdy2(i,j)+d2uwdz2(i,j));
                            nu*(d2uvdx2(i,j)+d2uvdy2(i,j)+d2uvdz2(i,j)) ...
                            nu*(d2vvdx2(i,j)+d2vvdy2(i,j)+d2vvdz2(i,j)) ...
                            nu*(d2vwdx2(i,j)+d2vwdy2(i,j)+d2vwdz2(i,j));
                            nu*(d2uwdx2(i,j)+d2uwdy2(i,j)+d2uwdz2(i,j)) ...
                            nu*(d2vwdx2(i,j)+d2vwdy2(i,j)+d2vwdz2(i,j)) ...
                            nu*(d2wwdx2(i,j)+d2wwdy2(i,j)+d2wwdz2(i,j))];
    end
end

VDxx=squeeze(squeeze(VD_tensor(1,1,:,:)));
VDyy=squeeze(squeeze(VD_tensor(2,2,:,:)));
VDzz=squeeze(squeeze(VD_tensor(3,3,:,:)));
VDxy=squeeze(squeeze(VD_tensor(1,2,:,:)));
VDxz=squeeze(squeeze(VD_tensor(1,3,:,:)));
VDyz=squeeze(squeeze(VD_tensor(2,3,:,:)));

%Velocity-pressure-gradient tensor. Tensor of Rank 2.
Pi_tensor(1:3,1:3,1:nx,1:ny)=0;

%Definition of the velocity-pressure-gradient tensor assuming 
%fully-developed flow, i.e., d()/dz=0, as a function of x and y.
%Piij=-1/rho*(<ui*dp/dxj>+<uj*dp/dxi>)
%Pi11=-2/rho*<u*dp/dx>
%Pi22=-2/rho*<v*dp/dy>
%Pi33=-2/rho*<w*dp/dz>
%Pi12=-1/rho*(<u*dp/dy>+<v*dp/dx>)
%Pi13=-1/rho*(<u*dp/dz>+<w*dp/dx>)
%Pi23=-1/rho*(<v*dp/dz>+<w*dp/dy>)

%Now, since we don't compute <ui*dp/dxj>, we use the chain rule to express
%these terms as a function of velocity gradients.
%<ui*dp/dxj>=d(<p*ui>)/dxj-<p*dui/dxj>

%We define the pressure transport and pressure strain tensors.

%Pressure transport tensor. Tensor of Rank 2.
PT_tensor(1:3,1:3,1:nx,1:ny)=0;

%Definition of the pressure transport tensor assuming 
%fully-developed flow, i.e., d()/dz=0, as a function of x and y.
%PTij=-1/rho*(<d(p*ui)/dxj>+<d(p*uj)/dxi>)
%PT11=-2/rho*<d(p*u)/dx>
%PT22=-2/rho*<d(p*v)/dy>
%PT33=-2/rho*<d(p*w)/dz>
%PT12=-1/rho*(<d(p*u)/dy>+<d(p*v)/dx>)
%PT13=-1/rho*(<d(p*u)/dz>+<d(p*w)/dx>)
%PT23=-1/rho*(<d(p*v)/dz>+<d(p*w)/dy>)

dpudx=D63-P.*dUdx-U.*D7;
dpvdx=D65-P.*dVdx-V.*D7;
dpwdx=D67-P.*dWdx-W.*D7;

dpudy=D64-P.*dUdy-U.*D8;
dpvdy=D66-P.*dVdy-V.*D8;
dpwdy=D68-P.*dWdy-W.*D8;

dpudz(1:nx,1:ny)=0;
dpvdz(1:nx,1:ny)=0;
dpwdz(1:nx,1:ny)=0;

for j=1:ny
    for i=1:nx
        PT_tensor(:,:,i,j)=[(dpudx(i,j)) ... 
                            (dpudy(i,j)+dpvdx(i,j)) ...
                            (dpudz(i,j)+dpwdx(i,j));
                            (dpudy(i,j)+dpvdx(i,j)) ...
                            (dpvdy(i,j)) ...
                            (dpvdz(i,j)+dpwdy(i,j));
                            (dpudz(i,j)+dpwdx(i,j)) ...
                            (dpvdz(i,j)+dpwdy(i,j)) ...
                            (dpwdz(i,j))];
    end
end

%Pressure strain tensor. Tensor of Rank 2.
PS_tensor(1:3,1:3,1:nx,1:ny)=0;

%Definition of the pressure strain tensor assuming 
%fully-developed flow, i.e., d()/dz=0, as a function of x and y.
%PSij=-1/rho*(<p*dui/dxj>+<p*duj/dxi>)
%PS11=-2/rho*<p*du/dx>
%PS22=-2/rho*<p*dv/dy>
%PS33=-2/rho*<p*dw/dz>
%PS12=-1/rho*(<p*du/dy>+<p*dv/dx>)
%PS13=-1/rho*(<p*du/dz>+<p*dw/dx>)
%PS23=-1/rho*(<p*dv/dz>+<p*dw/dy>)

%<pdudx>       % F15
%<pdudy>       % F16
%<pdudz>       % F17

%<pdvdx>       % F18
%<pdvdy>       % F19
%<pdvdz>       % F20

%<pdwdx>       % F21
%<pdwdy>       % F22
%<pdwdz>       % F23

pdudx=F15-P.*dUdx;
pdudy=F16-P.*dUdy;
pdudz=F17-P.*dUdz;

pdvdx=F18-P.*dVdx;
pdvdy=F19-P.*dVdy;
pdvdz=F20-P.*dVdz;

pdwdx=F21-P.*dWdx;
pdwdy=F22-P.*dWdy;
pdwdz=F23-P.*dWdz;

for j=1:ny
    for i=1:nx
        PS_tensor(:,:,i,j)=[(pdudx(i,j)) ... 
                            (pdudy(i,j)+pdvdx(i,j)) ...
                            (pdudz(i,j)+pdwdx(i,j));
                            (pdudy(i,j)+pdvdx(i,j)) ...
                            (pdvdy(i,j)) ...
                            (pdvdz(i,j)+pdwdy(i,j));
                            (pdudz(i,j)+pdwdx(i,j)) ...
                            (pdvdz(i,j)+pdwdy(i,j)) ...
                            (pdwdz(i,j))];
    end
end

%Construct velocity-pressure-gradient-tensor
for j=1:ny
    for i=1:nx
        Pi_tensor(:,:,i,j)=[-2/rho*(PT_tensor(1,1,i,j)-PS_tensor(1,1,i,j)) ...
                            -1/rho*(PT_tensor(1,2,i,j)-PS_tensor(1,2,i,j)) ...
                            -1/rho*(PT_tensor(1,3,i,j)-PS_tensor(1,3,i,j));
                            -1/rho*(PT_tensor(2,1,i,j)-PS_tensor(2,1,i,j)) ...
                            -2/rho*(PT_tensor(2,2,i,j)-PS_tensor(2,2,i,j)) ...
                            -1/rho*(PT_tensor(2,3,i,j)-PS_tensor(2,3,i,j));
                            -1/rho*(PT_tensor(3,1,i,j)-PS_tensor(3,1,i,j)) ...
                            -1/rho*(PT_tensor(3,2,i,j)-PS_tensor(3,2,i,j)) ...
                            -2/rho*(PT_tensor(3,3,i,j)-PS_tensor(3,3,i,j))];
    end
end

Pixx=squeeze(squeeze(Pi_tensor(1,1,:,:)));
Piyy=squeeze(squeeze(Pi_tensor(2,2,:,:)));
Pizz=squeeze(squeeze(Pi_tensor(3,3,:,:)));
Pixy=squeeze(squeeze(Pi_tensor(1,2,:,:)));
Pixz=squeeze(squeeze(Pi_tensor(1,3,:,:)));
Piyz=squeeze(squeeze(Pi_tensor(2,3,:,:)));


%Budget for each component of the Reynolds stress tensor 
%Without mean convection
Sxx=Pxx+Dxx+Txx+VDxx+Pixx;
Syy=Pyy+Dyy+Tyy+VDyy+Piyy;
Szz=Pzz+Dzz+Tzz+VDzz+Pizz;
Sxy=Pxy+Dxy+Txy+VDxy+Pixy;
Sxz=Pxz+Dxz+Txz+VDxz+Pixz;
Syz=Pyz+Dyz+Tyz+VDyz+Piyz;

%With mean convection
Scxx=Pxx+Dxx+Txx+VDxx+Pixx-Cxx;
Scyy=Pyy+Dyy+Tyy+VDyy+Piyy-Cyy;
Sczz=Pzz+Dzz+Tzz+VDzz+Pizz-Czz;
Scxy=Pxy+Dxy+Txy+VDxy+Pixy-Cxy;
Scxz=Pxz+Dxz+Txz+VDxz+Pixz-Cxz;
Scyz=Pyz+Dyz+Tyz+VDyz+Piyz-Cyz;

%Plot terms

if (ifplt)

%Mean velocities. Tensors of Rank 1.
figure(1)
pcolor(X,Y,U)
shading interp;colorbar;axis equal tight;

figure(2)
pcolor(X,Y,V)
shading interp;colorbar;axis equal tight;

figure(3)
pcolor(X,Y,W)
shading interp;colorbar;axis equal tight;

%Reynolds stress tensor. Tensor of Rank 2.
figure(4)
pcolor(X,Y,uu)
shading interp;colorbar;axis equal tight;

figure(5)
pcolor(X,Y,vv)
shading interp;colorbar;axis equal tight;

figure(6)
pcolor(X,Y,ww)
shading interp;colorbar;axis equal tight;

figure(7)
pcolor(X,Y,uv)
shading interp;colorbar;axis equal tight;

figure(8)
pcolor(X,Y,uw)
shading interp;colorbar;axis equal tight;

figure(9)
pcolor(X,Y,vw)
shading interp;colorbar;axis equal tight;

%Mean, RMS, skewness and flatness of pressure
figure(10)
pcolor(X,Y,P)
shading interp;colorbar;axis equal tight;

figure(11)
pcolor(X,Y,prms)
shading interp;colorbar;axis equal tight;

figure(12)
pcolor(X,Y,pskew)
shading interp;colorbar;axis equal tight;

figure(13)
pcolor(X,Y,pflat)
shading interp;colorbar;axis equal tight;

%Skewness
figure(14)
pcolor(X,Y,uuu)
shading interp;colorbar;axis equal tight;

figure(15)
pcolor(X,Y,vvv)
shading interp;colorbar;axis equal tight;

figure(16)
pcolor(X,Y,www)
shading interp;colorbar;axis equal tight;

figure(17)
pcolor(X,Y,uuv)
shading interp;colorbar;axis equal tight;

figure(18)
pcolor(X,Y,uuw)
shading interp;colorbar;axis equal tight;

figure(19)
pcolor(X,Y,uvv)
shading interp;colorbar;axis equal tight;

figure(20)
pcolor(X,Y,vvw)
shading interp;colorbar;axis equal tight;

figure(21)
pcolor(X,Y,uww)
shading interp;colorbar;axis equal tight;

figure(22)
pcolor(X,Y,vww)
shading interp;colorbar;axis equal tight;

figure(23)
pcolor(X,Y,uvw)
shading interp;colorbar;axis equal tight;

figure(24)
pcolor(X,Y,uskew)
shading interp;colorbar;axis equal tight;

figure(25)
pcolor(X,Y,vskew)
shading interp;colorbar;axis equal tight;

figure(26)
pcolor(X,Y,wskew)
shading interp;colorbar;axis equal tight;

% Production tensor
figure(27)
pcolor(X,Y,Pxx)
shading interp;colorbar;axis equal tight;

figure(28)
pcolor(X,Y,Pyy)
shading interp;colorbar;axis equal tight;

figure(29)
pcolor(X,Y,Pzz)
shading interp;colorbar;axis equal tight;

figure(30)
pcolor(X,Y,Pxy)
shading interp;colorbar;axis equal tight;

figure(31)
pcolor(X,Y,Pyz)
shading interp;colorbar;axis equal tight;

figure(32)
pcolor(X,Y,Pxz)
shading interp;colorbar;axis equal tight;

% Dissipation tensor
figure(33)
pcolor(X,Y,Dxx)
shading interp;colorbar;axis equal tight;

figure(34)
pcolor(X,Y,Dyy)
shading interp;colorbar;axis equal tight;

figure(35)
pcolor(X,Y,Dzz)
shading interp;colorbar;axis equal tight;

figure(36)
pcolor(X,Y,Dxy)
shading interp;colorbar;axis equal tight;

figure(37)
pcolor(X,Y,Dyz)
shading interp;colorbar;axis equal tight;

figure(38)
pcolor(X,Y,Dxz)
shading interp;colorbar;axis equal tight;

% Turbulent transport tensor
figure(39)
pcolor(X,Y,Txx)
shading interp;colorbar;axis equal tight;

figure(40)
pcolor(X,Y,Tyy)
shading interp;colorbar;axis equal tight;

figure(41)
pcolor(X,Y,Tzz)
shading interp;colorbar;axis equal tight;

figure(42)
pcolor(X,Y,Txy)
shading interp;colorbar;axis equal tight;

figure(43)
pcolor(X,Y,Tyz)
shading interp;colorbar;axis equal tight;

figure(44)
pcolor(X,Y,Txz)
shading interp;colorbar;axis equal tight;

% Turbulent transport tensor
figure(45)
pcolor(X,Y,VDxx)
shading interp;colorbar;axis equal tight;

figure(46)
pcolor(X,Y,VDyy)
shading interp;colorbar;axis equal tight;

figure(47)
pcolor(X,Y,VDzz)
shading interp;colorbar;axis equal tight;

figure(48)
pcolor(X,Y,VDxy)
shading interp;colorbar;axis equal tight;

figure(49)
pcolor(X,Y,VDyz)
shading interp;colorbar;axis equal tight;

figure(50)
pcolor(X,Y,VDxz)
shading interp;colorbar;axis equal tight;

% Velocity-pressure-gradient tensor
figure(51)
pcolor(X,Y,Pixx)
shading interp;colorbar;axis equal tight;

figure(52)
pcolor(X,Y,Piyy)
shading interp;colorbar;axis equal tight;

figure(53)
pcolor(X,Y,Pizz)
shading interp;colorbar;axis equal tight;

figure(54)
pcolor(X,Y,Pixy)
shading interp;colorbar;axis equal tight;

figure(55)
pcolor(X,Y,Piyz)
shading interp;colorbar;axis equal tight;

figure(56)
pcolor(X,Y,Pixz)
shading interp;colorbar;axis equal tight;

% Mean convection tensor
figure(57)
pcolor(X,Y,Cxx)
shading interp;colorbar;axis equal tight;

figure(58)
pcolor(X,Y,Cyy)
shading interp;colorbar;axis equal tight;

figure(59)
pcolor(X,Y,Czz)
shading interp;colorbar;axis equal tight;

figure(60)
pcolor(X,Y,Cxy)
shading interp;colorbar;axis equal tight;

figure(61)
pcolor(X,Y,Cyz)
shading interp;colorbar;axis equal tight;

figure(62)
pcolor(X,Y,Cxz)
shading interp;colorbar;axis equal tight;

% Balance without mean convection
figure(63)
pcolor(X,Y,Sxx)
shading interp;colorbar;axis equal tight;

figure(64)
pcolor(X,Y,Syy)
shading interp;colorbar;axis equal tight;

figure(65)
pcolor(X,Y,Szz)
shading interp;colorbar;axis equal tight;

figure(66)
pcolor(X,Y,Sxy)
shading interp;colorbar;axis equal tight;

figure(67)
pcolor(X,Y,Syz)
shading interp;colorbar;axis equal tight;

figure(68)
pcolor(X,Y,Sxz)
shading interp;colorbar;axis equal tight;

% Balance with mean convection
figure(69)
pcolor(X,Y,Scxx)
shading interp;colorbar;axis equal tight;

figure(70)
pcolor(X,Y,Scyy)
shading interp;colorbar;axis equal tight;

figure(71)
pcolor(X,Y,Sczz)
shading interp;colorbar;axis equal tight;

figure(72)
pcolor(X,Y,Scxy)
shading interp;colorbar;axis equal tight;

figure(73)
pcolor(X,Y,Scyz)
shading interp;colorbar;axis equal tight;

figure(74)
pcolor(X,Y,Scxz)
shading interp;colorbar;axis equal tight;

end

save(fsave)
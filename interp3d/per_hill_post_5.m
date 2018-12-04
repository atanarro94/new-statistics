close all; 
clear all; 
clc
double precision;
format long;

%% Load statistics
load 'data_stats.mat'

%% Select planes to be analysed (xy)
plane = ['xy'; 'yz'];
tol = 0.0001;

for j = 1:length(plane)
    if strcmp(plane(j,:),'xy')
        position = [0.5 1.5 2.5 3.5];
    elseif strcmp(plane(j,:),'xz')
        position = [0.5 1.5];
    elseif strcmp(plane(j,:),'yz')
        position = [-0.5 -1.5];
    end
    for i = 1:length(position)
        if strcmp(plane(j,:),'xy')
            vec = unique(z);
            ind{j}(:,i) = find(abs(z - position(i)) - min(abs(z-position(i))) < tol);
        elseif strcmp(plane(j,:),'xz')
            vec = unique(y);
            ind{j}(:,i) = find(abs(y - position(i)) - min(abs(y-position(i))) < tol);
        elseif strcmp(plane(j,:),'yz')
            vec = unique(x);
            ind{j}(:,i) = find(abs(x - position(i)) - min(abs(x-position(i))) < tol);
        end
    end
end

%% Wall units
Uinf = 1;
rho = 1;
mu = nu*rho;
pl = 'xy';

for i = 1:length(position)
    if strcmp(pl,'xz')
        indw = find(x == 0 & y == y(ind(1,i)));
        tauw(:,i) = mu*(dWdx(indw)+dVdx(indw));  % Friction in the wall at x = 0
    elseif strcmp(pl,'yz')
        indw = find(y == 0 & x == x(ind(1,i)));
        tauw(:,i) = mu*(dWdy(indw)+dUdy(indw));  % Friction in the wall at y = 0
    elseif strcmp(pl,'xy')  %TODO: adapt to use the average ut
        indw = find(x == 0 | y == 0);
        tauw(:,i) = mu*(dWdx(indw)+dVdx(indw)) + mu*(dWdy(indw)+dUdy(indw));
    end
end
ut = sqrt(abs(tauw)/rho);

%% Post-process data
pos = 3;
pla = 1;

% Coordinates
x1 = x(ind{pla}(:,pos));
y1 = y(ind{pla}(:,pos));
z1 = z(ind{pla}(:,pos));

% Mean flow
U1 = U(ind{pla}(:,pos));
V1 = V(ind{pla}(:,pos));
W1 = W(ind{pla}(:,pos));

% Reynolds stresses
uu1 = uu(ind{pla}(:,pos));
vv1 = vv(ind{pla}(:,pos));
ww1 = ww(ind{pla}(:,pos));
uv1 = uv(ind{pla}(:,pos));
uw1 = uw(ind{pla}(:,pos));
vw1 = vw(ind{pla}(:,pos));

P1 = V1;   % value to be plotted

if strcmp(plane(pla,:),'xy')
    xs = reshape(x1,[SX,SY]);
    ys = reshape(y1,[SX,SY]);
    Ps = reshape(P1,[SX,SY]);
elseif strcmp(plane(pla,:),'xz')
    xs = reshape(z1,[SX,SZ]);
    ys = reshape(x1,[SX,SZ]);
    Ps = reshape(P1,[SX,SZ]);
elseif strcmp(plane(pla,:),'yz')
    xs = reshape(z1,[SY,SZ]);
    ys = reshape(y1,[SY,SZ]);
    Ps = reshape(P1,[SY,SZ]);
end

figure
contourf(xs,ys,Ps,100,'LineColor','none')

% figure
% semilogx(ys(:,10)*ut_x(find(xp(:,10)==x1(1)))/nu,Ps(:,10)/ut_x(find(xp(:,10)==x1(1))))




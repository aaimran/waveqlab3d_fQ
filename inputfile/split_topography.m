close all
clc
clear all
nx = (351-1)/2 + 1;
ny = 183; %314/2;
nz = (761-1)/2 + 1;

% load fault profile
L = load('/home/kduru/Downloads/wasatch_100_surface_flat_curved_deep_topo.txt');

ax = 0;  bx = 60;

 figure(11)
Xc0k = reshape(L(:, 1)-5.378003800000000, ny, nz);
Yc0k = reshape(L(:, 2)-0.5806683401649992, ny, nz);
Zc0k = reshape(L(:, 3), ny, nz);
% 
% 
mesh(Zc0k, Yc0k, Xc0k+30)
xlim([0, 60])
ylim([-2, 30])
zlim([14, 30])

% extract boundary curves for the first block
XR(1:nz) = Xc0k(1, :)+30; XL(1:nz) = 0;

h1x = (XR(1)-ax)/(nx-1);
hnx = (XR(end)-ax)/(nx-1);



XB(1:nx) = 0:h1x:XR(1); XT(1:nx) = 0:hnx:XR(end);

[xleft, xright] = stretchedgrids(nx);

XB(1:nx) = 0 + (XR(1)-0)*xleft; XT(1:nx) = 0 + (XR(end)-0)*xleft;

% generate a curvilinear mesh obeying the boundary curves
[X1] = Gen_Surf_Mesh(nx,nz, XL,XR,XB,XT);

figure(12)
mesh(X1)

%return

% extract boundary curves for the first block
XL(1:nz) = Xc0k(1, :)+30; XR(1:nz) = 60;

h1x = (bx-XL(1))/(nx-1);
hnx = (bx-XL(end))/(nx-1);

XB(1:nx) = XL(1):h1x:bx; XT(1:nx) = XL(end):hnx:bx;

XB(1:nx) = XL(1) + (bx-XL(1))*xright; XT(1:nx) = XL(end) + (bx-XL(end))*xright;


% generate a curvilinear mesh obeying the boundary curves
[X2] = Gen_Surf_Mesh(nx,nz, XL,XR,XB,XT);

figure(13)
mesh(X2)

%return
% merge both blocks for a single mesh
X = [X1 X2(:, 2:end)];

figure(14)
mesh(X)

%return
% load topography data for the first block
K = load('topo_block1_dec.txt');
nx = (351-1)/2 + 1;
ny =  314/2;
nz = (761-1)/2 + 1;

X1k = reshape(K(:, 1), nz, nx);
Y1k = reshape(K(:, 2), nz, nx);
Z1k = reshape(K(:, 3), nz, nx);

figure(15)
mesh(X1k, Z1k, Y1k)
axis image
set(gca,'zdir','reverse')
xlabel('fault normal')
ylabel('fault parallel')
zlabel('along dip')

% size(Xk)
% max(max(abs(Yk)))

figure(20)

plot(Z1k(:, end), Y1k(:, end))
%hold on 


% load topography data for the second block
K = load('topo_block2_dec.txt');

nx = (351-1)/2 + 1;
ny = 183; %314/2;
nz = (761-1)/2 + 1;

X2k = reshape(K(:, 1), nz, nx);
Y2k = reshape(K(:, 2), nz, nx);
Z2k = reshape(K(:, 3), nz, nx);
figure(16)

mesh(X2k, Z2k, Y2k)
axis image
set(gca,'zdir','reverse')
xlabel('fault normal')
ylabel('fault parallel')
zlabel('along dip')

% combine to a single block
X0 = [X1k X2k(:, 2:end)];
Y0 = [Y1k Y2k(:, 2:end)];
Z0 = [Z1k Z2k(:, 2:end)];

figure(17)
mesh(X0, Z0, Y0)

axis image
set(gca,'zdir','reverse')
xlabel('fault normal')
ylabel('fault parallel')
zlabel('along dip')

[mz, my] = size(X0);

Y = zeros(mz, my);
Z = zeros(mz, my);

% write data to a file
write_topography(X0,Y0,Z0,'curv_topo', mz, my)

XX = load('curv_topo.dat');

XX0 = reshape(XX(:, 1), mz, my);
YY0 = reshape(XX(:, 2), mz, my);
ZZ0 = reshape(XX(:, 3), mz, my);


% Contruct high order interpolation such that topography data 
% matches fault trace on the surface
%Xc0k-5.378003800000000

%return
for i = 1:mz
    
    for j = 1:my
        
        x = X(i, j);
        z = Z0(i, j);   

        [y] = interpolate_surface(Y0, Z0, X0, z, x, mz, my);  
        
        Z(i, j) = z;
        Y(i, j) = y;
        
    end
end

% split the topography data along the fault trace into two equal files
X1 = X(:, 1:nx);
Y1 = Y(:, 1:nx);
Z1 = Z(:, 1:nx);

X2 = X(:, nx:end);
Y2 = Y(:, nx:end);
Z2 = Z(:, nx:end);

figure(4)
mesh(X1, Z1, Y1)
set(gca,'zdir','reverse')
xlabel('fault normal')
ylabel('fault parallel')
zlabel('along dip')
axis image
hold on

%figure(5)
mesh(X2, Z2, Y2)
set(gca,'zdir','reverse')
xlabel('fault normal')
ylabel('fault parallel')
zlabel('along dip')
axis image

hold off


% write data to a file
write_topography(X1,Y1,Z1,'block1_curv_topo_stretched', nz,nx)
write_topography(X2,Y2,Z2,'block2_curv_topo_stretched', nz,nx)


K1 = load('block1_curv_topo_stretched.dat');


X1k = reshape(K1(:, 1), nz, nx);
Y1k = reshape(K1(:, 2), nz, nx);
Z1k = reshape(K1(:, 3), nz, nx);

K2 = load('block2_curv_topo_stretched.dat');


X2k = reshape(K2(:, 1), nz, nx);
Y2k = reshape(K2(:, 2), nz, nx);
Z2k = reshape(K2(:, 3), nz, nx);

figure(41)
mesh(X1k, Z1k, Y1k)
set(gca,'zdir','reverse')
xlabel('fault normal')
ylabel('fault parallel')
zlabel('along dip')
axis image
hold on

%figure(5)
mesh(X2k, Z2k, Y2k)
set(gca,'zdir','reverse')
xlabel('fault normal')
ylabel('fault parallel')
zlabel('along dip')
axis image

hold off




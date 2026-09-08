clc
close all
clear all

% must match by in the inputfile
nz = 201;
nx = 101;
ny = 151;
by = 15;     

n0y = 171;    % provided by data

% load topography data for the first block
D = load('~/Downloads/block1_topo.dat');

X1k = reshape(D(:, 1), nz, nx);
Y1k = reshape(D(:, 2), nz, nx);
Z1k = reshape(D(:, 3), nz, nx);

figure(4)
mesh(X1k, Z1k, Y1k)
set(gca,'zdir','reverse')
xlabel('fault normal')
ylabel('fault parallel')
zlabel('along dip')
%axis image
hold on


% load topography data for the second block
F = load('~/Downloads/block2_topo.dat');

X2k = reshape(F(:, 1), nz, nx);
Y2k = reshape(F(:, 2), nz, nx);
Z2k = reshape(F(:, 3), nz, nx);

%figure(5)
mesh(X2k, Z2k, Y2k)
set(gca,'zdir','reverse')
xlabel('fault normal')
ylabel('fault parallel')
zlabel('along dip')
%axis image

hold off

% load fault surface profile
C = load('~/Downloads/fault_topo.dat');

Xc0k = reshape(C(:, 1), n0y, nz);
Yc0k = reshape(C(:, 2), n0y, nz);
Zc0k = reshape(C(:, 3), n0y, nz);


figure(3)
mesh(Zc0k, Yc0k, Xc0k)
axis image
set(gca,'zdir','reverse')
xlabel('fault normal')
ylabel('fault parallel')
zlabel('along dip')

% need generate interpolated fault surface constrained against topography


YT(1:nz) = Y1k(:, end);
YB(1:nz) = by; 

r = linspace(0, 1, ny);

hy1 = (YB(1) - YT(1))/(ny-1);
hyn = (YB(end) - YT(end))/(ny-1);

YL(1:ny) = YT(1):hy1:YB(1);
YR(1:ny) = YT(end):hyn:YB(end);
YL(1:ny) = YT(1) + (YB(1)-YT(1))*r;
YR(1:ny) = YT(end) + (YB(end)-YT(end))*r;

ZT(1:nz) = Z1k(:, end);
ZB(1:nz) = Zc0k(end, :);
ZL(1:ny) = ZT(1);
ZR(1:ny) = ZT(end);
    

[Y1] = Gen_Surf_Mesh(nz,ny, YL,YR,YT,YB);
[Z1] = Gen_Surf_Mesh(nz,ny, ZL,ZR,ZT,ZB);


X = zeros(ny, nz);
Y = zeros(ny, nz);
Z = zeros(ny, nz);


for i = 1:ny
    
    for j = 1:nz
        
        y = Y1(i, j);
        z = Z1(i, j);   

        [x] = interpolate_surface(Xc0k, Yc0k, Zc0k, y, z, 171, nz);  
        
        X(i, j) = x;
        
        Y(i, j) = y;
        
        Z(i, j) = z;
        
    end
end

drawnow
figure(43)
mesh(Z, Y, X)

write_fault_topography(X',Y',Z','Gaussian_100m', nz, ny)



L = load('Gaussian_100m.dat');



Xck = reshape(L(:, 1), ny, nz);
Yck = reshape(L(:, 2), ny, nz);
Zck = reshape(L(:, 3), ny, nz);


drawnow
figure(44)
mesh(Zck, Yck, Xck)





% test_least_squares_2d_iter_LT.m

close all;
clear;
clc;

x = 0 : 1 : 10;
y = 0 : 1 : 10;
[gx, gy] = meshgrid(x, y);

f = sin(gx).*cos(gy);

niter = 3;
Lf = f; h = 1;
for iter = 1 : niter
    Lf = least_squares_2d(Lf, h);
    h = 1/2^iter;
end

x1 = 0 : h : 10;
y1 = 0 : h : 10;
[gx1, gy1] = meshgrid(x1, y1);

figure; 
surf(gx, gy, f);

figure;
surf(gx1, gy1, Lf);




function Lf = least_squares_2d(f, h)

% Find kernels.
x = (0 : 1 : 3)*h;
y = (0 : 1 : 3)*h;
[gx, gy] = meshgrid(x, y);

U = [gx(:).^0, gx(:), gy(:), gx(:).^2, gx(:).*gy(:), gy(:).^2, gx(:).^3, gx(:).^2.*gy(:), ...
     gx(:).*gy(:).^2, gy(:).^3];

f_pad = [f(:,1), f, f(:,end), f(:,end)];
f_pad = [f_pad(1,:); f_pad; f_pad(end,:); f_pad(end,:)];

% vertex points
xb = 1*h; yb = 1*h;
W = diag(exp(-((gx(:)-xb).^2 + (gy(:)-yb).^2)));
u = [xb.^0, xb, yb, xb^2, xb*yb, yb^2, xb^3, xb^2*yb, xb*yb^2, yb^3]';
L = W*U*inv(U'*W*U)*u;
L = reshape(L, [4,4]);
Lf_vert = filter2(L, f_pad, 'valid');

% horizontal edge points
xb = 3/2*h; yb = 1*h;
W = diag(exp(-((gx(:)-xb).^2 + (gy(:)-yb).^2)));
u = [xb.^0, xb, yb, xb^2, xb*yb, yb^2, xb^3, xb^2*yb, xb*yb^2, yb^3]';
L = W*U*inv(U'*W*U)*u;
L = reshape(L, [4,4]);
Lf_hedge = filter2(L, f_pad, 'valid');

% vertical edge points
xb = 1*h; yb = 3/2*h;
W = diag(exp(-((gx(:)-xb).^2 + (gy(:)-yb).^2)));
u = [xb.^0, xb, yb, xb^2, xb*yb, yb^2, xb^3, xb^2*yb, xb*yb^2, yb^3]';
L = W*U*inv(U'*W*U)*u;
L = reshape(L, [4,4]);
Lf_vedge = filter2(L, f_pad, 'valid');

% face points
xb = 3/2*h; yb = 3/2*h;
W = diag(exp(-((gx(:)-xb).^2 + (gy(:)-yb).^2)));
u = [xb.^0, xb, yb, xb^2, xb*yb, yb^2, xb^3, xb^2*yb, xb*yb^2, yb^3]';
L = W*U*inv(U'*W*U)*u;
L = reshape(L, [4,4]);
Lf_face = filter2(L, f_pad, 'valid');

[R, C] = size(f);
Lf = zeros(2*R-1, 2*C-1);
Lf(1:2:end, 1:2:end) = Lf_vert;
Lf(1:2:end, 2:2:end) = Lf_hedge(:,1:end-1);
Lf(2:2:end, 1:2:end) = Lf_vedge(1:end-1,:);
Lf(2:2:end, 2:2:end) = Lf_face(1:end-1,1:end-1);

end







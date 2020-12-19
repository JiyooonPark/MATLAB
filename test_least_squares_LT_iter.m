
% test_least_squares_iter_LT.m

% Iterative least squares method

close all;
clear;
clc;

x = 0 : 1 : 10;
f = sin(x);

figure;
plot(x, f, 'o-');

niter = 3;
Lf = f; h = 1;
for iter = 1 : niter
    Lf = least_squares(Lf, h);
    h = 1/2^iter;
end


x1 = 0 : h : 10;
figure; hold on;
plot(x, f, 'o-');
plot(x1, Lf, 'o-');


%%
function Lf = least_squares(f, h)

    x = (0 : 1 : 5)*h;
    U = [x(:).^0, x(:), x(:).^2, x(:).^3];

    % old
    xb = 2*h;
    W = diag(exp(-(x(:)-xb).^2));
    u = [xb.^0, xb, xb^2, xb^3]';
    L = W*U*inv(U'*W*U)*u
    f_pad = [f(1) f(1) f, f(end) f(end) f(end)];    
    Lf_old = conv(f_pad, L(end:-1:1)', 'valid');

    % new
    xb = 5/2*h;
    W = diag(exp(-(x(:)-xb).^2));
    u = [xb.^0, xb, xb^2, xb^3]';
    L = W*U*inv(U'*W*U)*u       
    f_pad = [f(1) f(1) f, f(end) f(end)];
    Lf_new = conv(f_pad, L(end:-1:1)', 'valid');
    
    N = length(f);
    Lf = zeros(1, 2*N-1);
    Lf(1:2:end) = Lf_old;
    Lf(2:2:end) = Lf_new;

end









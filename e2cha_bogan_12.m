close all;
clc;
clear;

x0 = 0;
x1=1;
x2=4;


x3 = (f(x0) * (x1^2 - x2^2) + f(x1) * (x2^2 - x0^2) ...
    + f(x2) * (x0^2 - x1^2))/ (2*f(x0) * (x1-x2) + ...
    2* f(x1)* (x2-x0) + 2* f(x2) * (x0-x1))
f(x3)
figure; hold on;
x=[x0 x1 x2 x3];
y=[f(x0), f(x1), f(x2), f(x3)];
plot(x, y, 'o:');
axis equal;


function xr = f(x)
    xr = 2*sin(x)-x^2/10;
end

% test_SVD.m

close all;
clear;
clc;
% 
% A = sym( [1, 0;
%           0, 1;
%           1, 1] )
%       % 숫자가 아니라 기호로 인식을 해라 
%       
% [V, D] = eig(A'*A)
% 
% V = [V(:, 2), (-1)*V(:, 1)]
% D = diag([D(2,2), D(1,1)])
% 
% V(:, 1) = V(:, 1) / norm(V(:, 1))
% V(:, 2) = V(:, 2) / norm(V(:, 2))
% 
% u1 = simplify( A*V(:,1) / sqrt(D(1,1)) )
% u2 = simplify( A*V(:,2) / sqrt(D(2,2)) )
% 
% [U, E] = eig(A*A')
% u3 = (-1)*U(:, 1)
% u3 = u3 / norm(u3) % 노멀라이즈 한 것이다. 
% 
% U = [u1, u2, u3]
% S = [sqrt(D(1,1)), 0;
%      0, sqrt(D(2,2));
%      0, 0] %sigma
%  
% A1 = simplify( U*S*V' )

A = [1, 0;
     0, 1;
     1, 1]
[U, S, V] = svd(A)




 
 
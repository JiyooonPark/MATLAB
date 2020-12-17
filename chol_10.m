close all;
clc;
clear;

A=[6 15 55;
    15 55 225;
    55 225 979;];
R = chol(A)

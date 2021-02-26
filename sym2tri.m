clc;
clear all;
close all;

n = input('Enter size of matrix - ');
% Generate a symmetric matrix
A = rand(n, n);
A = (A + transpose(A))/2;
% hessenberg is a function used to convert matrix to a hessenberg
% matrix
H = hessenberg(A, n);
% disp(H);

clc; clear all; close all;

% Δημιουργία Λευκών θορύβων για να έχουμε 
% τα ίδια αποτελέσματα στην
% προσομοίωση
q = 2.5; r = 1.3;
N = 500;
w = sqrt(q)*randn(N,1);
v = sqrt(r)*randn(N,1);
save('w.mat','w');
save('v.mat','v');
clc; clear all; close all;

% Section A
%%
% Question 1
syms s;
%������ ��� �������� ��� ���������� ���������
S = (s-20)/((s+500)*(s+200)*(s+100)*(s+2)*(s+1));
S = expand(S);
S = simplify(S);
pretty(S);
%%
% Question 2
%���������� ��� ���������� ���������
num = [1 -20];
den = [1 803 172402 10511600 30340000 20000000];
Gs = tf(num,den);
%%
% ���������� Observable Cannonical Form
% ���������� ���� �� ����� ������������ �� �������
n = [0 0 0 1 -20];
d = [803 172402 10511600 30340000 20000000];
A = [0 0 0 0 -d(1,5);
     1 0 0 0 -d(1,4);
     0 1 0 0 -d(1,3);
     0 0 1 0 -d(1,2);
     0 0 0 1 -d(1,1)];
B = [n(1,5); n(1,4); n(1,3); n(1,2); n(1,1)];
C = [0 0 0 0 1];
D = 0;
%%
% ��������� �� Jordan Cannonical Form;
[V,J] = jordan(A);
Bj = inv(V)*B;
Cj = C*V;
J = diag([-500 -200 -100 -2 -1]);
BJ = [Bj(2,1) ; Bj(3,1) ; Bj(1,1) ; Bj(4,1) ; Bj(5,1)];
sys = ss(J,BJ,Cj,D);
%%
%Question 3 - 4
% ������� �������
G = [1 ; 0 ; 0 ; 1 ; 1];
W = [1 ; 1 ; 1];
% ������� ����������� ���������
H = [ 1 0 0 0 0 ;
      0 0 1 0 0 ;
      0 0 0 0 1]; 
% ���� ����� ����������
sysA = ss(J,[Bj G 0*Bj] , [Cj ; H] , [D zeros(1,2); zeros(3,2) W]);
%%
% ���������� ����������� ���������
[n1,d1] = ss2tf(sysA.A,sysA.B,sysA.C,sysA.D,1);
[n2,d2] = ss2tf(sysA.A,sysA.B,sysA.C,sysA.D,2);
[n3,d3] = ss2tf(sysA.A,sysA.B,sysA.C,sysA.D,3);
tf11 = tf(n1(1,:),d1); tf12 = tf(n1(2,:),d1); 
tf13 = tf(n1(3,:),d1); tf14 = tf(n1(4,:),d1);
tf21 = tf(n2(1,:),d2); tf22 = tf(n2(2,:),d2);
tf23 = tf(n2(3,:),d2); tf24 = tf(n2(4,:),d2);
tf31 = tf(n3(1,:),d3); tf32 = tf(n3(2,:),d3); 
tf33 = tf(n3(3,:),d3); tf34 = tf(n3(4,:),d3);
% ������ Bandwidth ��� ��� ����������� ���������
BW = [bandwidth(tf11) bandwidth(tf12) bandwidth(tf13) bandwidth(tf14);
      bandwidth(tf21) bandwidth(tf22) bandwidth(tf23) bandwidth(tf24);
      bandwidth(tf31) bandwidth(tf32) bandwidth(tf33) bandwidth(tf34) ];
% ������� ��������� Bandwidth
wbw = 0.8391;
% ������� �������� ��������������
Tsmax = (2*pi)/(35*wbw);
%%
% ������� �������� �������������� 
Ts = 0.1;
% �������������� ����������
sysAd = c2d(sysA , Ts , 'zoh');
% ���������������� 
if rank(obsv(sysAd.A,sysAd.C)) == 5 
     fprintf('The discrete system is Observable \n');
end
% ��������� ����������
eig(sysAd);
%%
% Question 5
% ������������ �������
q = 2.5; r = 1.3;
% ������� ������������ 
Q = q*eye(5); R = r*eye(3);
% ������������ ����������
x0 = [0 ; 0 ; 0 ; 0 ; 0];
P0 = eye(5);
% ������� �������� ��������� ��� ������������
% ��� ����������
t = (0:Ts:50-Ts);
t = t';
% ������� ��� ���������� 
u = sin(2*t);
% �� ���������� ��� � ������� ��� ���������
load w.mat;
load v.mat;
% ����������� ��� ����������
[out,~,statesA] = lsim(sysAd,[u  w  v],t,x0);
out = out' ;
% �������� ���������
z = [ out(2,:) ; out(3,:) ; out(4,:) ];
%%
% Question 6
% �������� ������������� ����� ��� ������ ���
% ������� �� ������������
innovationA = zeros(3,length(t));
x_eA = zeros(5,length(t));
% ������� � ��� ������� ����������
Bd = sysAd.B(:,1);
% ������� ���������� ��� ��� ��������� ��� ��� ��� 
% ���������������� ���������
x = x0;
P = P0;
% ��������� ���������� ��� ������� Kalman
for i=1:length(t)
    %��������
    %x[k](a priori) = Ax[k-1] + Bu[k-1]
    x = sysAd.A*x+Bd*u(i,:); 
    % P[k](a priori) = A*P[k-1]*A' + Q
    P = A*P*A'+Q;  
    %��������
    % K = P[k](a priori)*H'/(H*P[k](a priori)*H')
    K = P*H'/(H*P*H'+R);
    % ����������� �����������
    innovationA(:,i) = z(:,i) - H*x ;
    % x[k](a posteriori) = x[k](a priori) + K*(z[k]-H*x[k](a priori))
    x = x + K*innovationA(:,i);
    % P[k](a posteriori) = (I-K*H)*P(a priori)
    P = (eye(5)-K*H)*P;
    % �������� ��������� ��� �����������
    x_eA(:,i) = x; 
end
%%
% �������� ����� ��� ���������� �������� 
% ���������
error = zeros(5,length(t));
N=length(t);
% ������� ��� ��������� ��� ������� Kalman
for j = 1:5
    % ������ ��������� ������ ��� ����������� 
    % ��� ��� ��������� �����
    error(j,:) = statesA(:,j) - (x_eA(j,:))';
end
%%
% �������� ���������� ���������
plot(t,error(1,:));
grid minor;
title('Estimation Error of x1');
xlabel('Time(s)');
pause();
close;

plot(t,error(2,:));
grid minor;
title('Estimation Error of x2');
xlabel('Time(s)');
pause();
close;

plot(t,error(3,:));
grid minor;
title('Estimation Error of x3');
xlabel('Time(s)');
pause();
close;

plot(t,error(4,:));
grid minor;
title('Estimation Error of x4');
xlabel('Time(s)');
pause();
close;

plot(t,error(5,:));
grid minor;
title('Estimation Error of x5');
xlabel('Time(s)');
pause();
close;
%%
% ����� ���� �������� ������� Kalman
% ��� ���� ���������
p = zeros(5,1);
for j = 1:5
    % �������������� �� ���� ��� ������� ��������
    error(j,:) = error(j,:)/max(abs(error(j,:)));
    % ������ ������ ��� ����������� ���������
    p(j) = norm(error(j,:));
end
% ����� RMSE 
performance = 1- sqrt(sum(p)/N)
%%
% ������������ �����������
% �������������� ��� ����������� ���� 
% �� ����� ����� ��� [-1 ,1] ��� �� ����� �������� 
% � ���� ���� 
innovationA1_n = innovationA(1,:)/max(abs(innovationA(1,:)));
innovationA2_n = innovationA(2,:)/max(abs(innovationA(2,:)));
innovationA3_n = innovationA(3,:)/max(abs(innovationA(3,:)));
% ����� ����� ���������� 
means_in = [ sum(innovationA1_n)/N ; 
             sum(innovationA2_n)/N ;
             sum(innovationA3_n)/N ];
% ����� ����� ������ ������� ��� ��������
mean_w = sum(w)/length(w);
mean_v = sum(v)/length(v);
mean_noise = [ mean_w mean_v];
%%
% ���������� ������������� ��� �� ���������� �� 
% ����������� Gaussian ��������
histogram(innovationA1_n)
title('Innovation of x1');
pause();
close;
histogram(innovationA2_n)
title('Innovation of x3');
pause();
close;
histogram(innovationA3_n)
title('Innovation of x5');
pause();
close;
%%
save('statesA.mat','statesA');
save('x_eA.mat','x_eA');


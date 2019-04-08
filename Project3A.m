clc; clear all; close all;

% Section A
%%
% Question 1
syms s;
%Εύρεση της εξίσωσης της συνάρτησης μεταφοράς
S = (s-20)/((s+500)*(s+200)*(s+100)*(s+2)*(s+1));
S = expand(S);
S = simplify(S);
pretty(S);
%%
% Question 2
%Δημιουργία της συνάρτησης μεταφοράς
num = [1 -20];
den = [1 803 172402 10511600 30340000 20000000];
Gs = tf(num,den);
%%
% Δημιουργια Observable Cannonical Form
% συστήματος ώστε να είναι παρατηρήσιμο το σύστημα
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
% Μετατροπή σε Jordan Cannonical Form;
[V,J] = jordan(A);
Bj = inv(V)*B;
Cj = C*V;
J = diag([-500 -200 -100 -2 -1]);
BJ = [Bj(2,1) ; Bj(3,1) ; Bj(1,1) ; Bj(4,1) ; Bj(5,1)];
sys = ss(J,BJ,Cj,D);
%%
%Question 3 - 4
% Πίνακας Θορυβων
G = [1 ; 0 ; 0 ; 1 ; 1];
W = [1 ; 1 ; 1];
% Πίνακας ευαισθησίας μετρήσεων
H = [ 1 0 0 0 0 ;
      0 0 1 0 0 ;
      0 0 0 0 1]; 
% Νέος χώρος κατάστασης
sysA = ss(J,[Bj G 0*Bj] , [Cj ; H] , [D zeros(1,2); zeros(3,2) W]);
%%
% Σχεδιασμός συναρτήσεων μεταφοράς
[n1,d1] = ss2tf(sysA.A,sysA.B,sysA.C,sysA.D,1);
[n2,d2] = ss2tf(sysA.A,sysA.B,sysA.C,sysA.D,2);
[n3,d3] = ss2tf(sysA.A,sysA.B,sysA.C,sysA.D,3);
tf11 = tf(n1(1,:),d1); tf12 = tf(n1(2,:),d1); 
tf13 = tf(n1(3,:),d1); tf14 = tf(n1(4,:),d1);
tf21 = tf(n2(1,:),d2); tf22 = tf(n2(2,:),d2);
tf23 = tf(n2(3,:),d2); tf24 = tf(n2(4,:),d2);
tf31 = tf(n3(1,:),d3); tf32 = tf(n3(2,:),d3); 
tf33 = tf(n3(3,:),d3); tf34 = tf(n3(4,:),d3);
% Εύρεση Bandwidth των δυο συναρτήσεων μεταφοράς
BW = [bandwidth(tf11) bandwidth(tf12) bandwidth(tf13) bandwidth(tf14);
      bandwidth(tf21) bandwidth(tf22) bandwidth(tf23) bandwidth(tf24);
      bandwidth(tf31) bandwidth(tf32) bandwidth(tf33) bandwidth(tf34) ];
% Επιλογή ελαχίστου Bandwidth
wbw = 0.8391;
% Μεγιστη περίοδος δειγματοληψίας
Tsmax = (2*pi)/(35*wbw);
%%
% Επιλογή περιόδου δειγματοληψίας 
Ts = 0.1;
% Διακριτοποίηση σύστηματος
sysAd = c2d(sysA , Ts , 'zoh');
% Παρατηρησιμότητα 
if rank(obsv(sysAd.A,sysAd.C)) == 5 
     fprintf('The discrete system is Observable \n');
end
% Ιδιοτιμές συστήματος
eig(sysAd);
%%
% Question 5
% Συνδιασπορές Θορύβων
q = 2.5; r = 1.3;
% Πίνακες συνδιασπορων 
Q = q*eye(5); R = r*eye(3);
% Αρχικοποίηση παραμέτρων
x0 = [0 ; 0 ; 0 ; 0 ; 0];
P0 = eye(5);
% Χρονικό διάστημα εκτέλεσης της προσομοίωσης
% του συστήματος
t = (0:Ts:50-Ts);
t = t';
% Είσοδος του συστήματος 
u = sin(2*t);
% Οι διαταραχές και ο θόρυβος των μετρήσεων
load w.mat;
load v.mat;
% Προσομοίωση του συστήματος
[out,~,statesA] = lsim(sysAd,[u  w  v],t,x0);
out = out' ;
% Διάνυσμα μετρήσεων
z = [ out(2,:) ; out(3,:) ; out(4,:) ];
%%
% Question 6
% Δέσμευση αποθηκευτικού χώρου για μεγέθη που
% θέλουμε να υπολογίσουμε
innovationA = zeros(3,length(t));
x_eA = zeros(5,length(t));
% Πίνακας Β του αρχικού συστήματος
Bd = sysAd.B(:,1);
% Αρχικές Εκτιμήσεις για την κατάσταση και για την 
% συμμεταβλητότητα σφάλματος
x = x0;
P = P0;
% Υλοποίηση Αλγορίθμου του φίλτρου Kalman
for i=1:length(t)
    %Πρόβλεψη
    %x[k](a priori) = Ax[k-1] + Bu[k-1]
    x = sysAd.A*x+Bd*u(i,:); 
    % P[k](a priori) = A*P[k-1]*A' + Q
    P = A*P*A'+Q;  
    %Διόρθωση
    % K = P[k](a priori)*H'/(H*P[k](a priori)*H')
    K = P*H'/(H*P*H'+R);
    % Υπολογισμός καινοτομιών
    innovationA(:,i) = z(:,i) - H*x ;
    % x[k](a posteriori) = x[k](a priori) + K*(z[k]-H*x[k](a priori))
    x = x + K*innovationA(:,i);
    % P[k](a posteriori) = (I-K*H)*P(a priori)
    P = (eye(5)-K*H)*P;
    % Διάνυσμα εκτίμησης των καταστάσεων
    x_eA(:,i) = x; 
end
%%
% Δέσμευση χώρου για μεγαλύτερη ταχύτητα 
% εκτέλεσης
error = zeros(5,length(t));
N=length(t);
% Απόδοση της εκτίμησης του φίλτρου Kalman
for j = 1:5
    % Εύρεση σφάλματος μεταξύ των καταστάσεων 
    % και της εκτίμησης αυτών
    error(j,:) = statesA(:,j) - (x_eA(j,:))';
end
%%
% Γραφικές Παρστάσεις σφαλμάτων
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
% Μέσος όρος απόδοσης φίλτρου Kalman
% για κάθε κατάσταση
p = zeros(5,1);
for j = 1:5
    % Κανονικοποίηση ως προς την μέγιστη απόκλιση
    error(j,:) = error(j,:)/max(abs(error(j,:)));
    % Εύρεση νόρμας των διανυσμάτων σφάλματος
    p(j) = norm(error(j,:));
end
% Τύπος RMSE 
performance = 1- sqrt(sum(p)/N)
%%
% Απεικονίσεις Καινοτομιών
% Κανονικοποίηση των καινοτομιών ώστε 
% να έχουν εύρος απο [-1 ,1] και να φανεί καλύτερα 
% η μέση τιμή 
innovationA1_n = innovationA(1,:)/max(abs(innovationA(1,:)));
innovationA2_n = innovationA(2,:)/max(abs(innovationA(2,:)));
innovationA3_n = innovationA(3,:)/max(abs(innovationA(3,:)));
% Μέσες τιμές καινοτομών 
means_in = [ sum(innovationA1_n)/N ; 
             sum(innovationA2_n)/N ;
             sum(innovationA3_n)/N ];
% Μέσες τιμές λευκών θορύβων για συγκριση
mean_w = sum(w)/length(w);
mean_v = sum(v)/length(v);
mean_noise = [ mean_w mean_v];
%%
% Απεικόνιση ιστογραμμάτων για να εξετάσουμε αν 
% ικανοποιούν Gaussian κατανομή
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


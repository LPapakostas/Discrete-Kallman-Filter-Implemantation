clc; close all; clear all;
% ξανααααααα
%Question 1
% Ορισμός του πλήρες συστήματος
A = diag([-500 -200 -100 -2 -1]);
B = [-1.744e-8 ; 1.861e-7 ; -3.092e-7 ; 2.277e-6 ; -2.136e-6];
% Νέος πίνακας ευαισθησίας των μετρήσεων
H = [0 0 0 1 0 ; 0 0 0 0 1];
% Νέοι πίνακες θορύβων
G = [0 ; 0 ; 0 ; 1 ; 1];
W = [1 ; 1 ];
% Νέος χώρος κατάστασης 
sysB = ss(A,[B G 0*B] , H , [zeros(2,1) zeros(2,1) W]);
%%
% Επιλογή περιόδου δειγματοληψίας 
Ts = 0.1;
% Διακριτοποίηση πλήρους συστήματος
sysBd = c2d(sysB,Ts,'zoh');
%%
% Κρατάμε μόνο το υποσύστημα με τις 2 τελευταίες καταστάσεις
A1 = [sysB.A(4,4) sysB.A(4,5);
      sysB.A(5,4) sysB.A(5,5)];
B1 = [sysB.B(4,1) ; sysB.B(5,1)];
G1 = [sysB.B(4,2) ; sysB.B(5,2)];
H1 = [sysB.C(1,4) sysB.C(1,5);
      sysB.C(2,4) sysB.C(2,5)];
W1 = sysB.D(:,3);
% Μειωμένο Σύστημα
sysm = ss(A1,[B1 G1 0*B1],H1,[0*W1 0*W1 W1]);
% Διακριτοποίηση μειωμένου συστήματος
sysmd = c2d(sysm , Ts , 'zoh');
% Παρατηρησιμότητα 
if rank(obsv(sysmd.A,sysmd.C)) == 2
     fprintf('The discrete subsystem is Observable \n');
end
%%
% Question 2
% Συνδιασπορές θορύβων
q1 = 5;  r = 1.3;
% Πίνακες Συνδιασπορών
Q1 = q1*eye(2);
R = r*eye(2);
%%
% Χρονικό διάστημα εκτέλεσης της προσομοίωσης
% του συστήματος
t = (0:Ts:50-Ts);
t = t';
% Είσοδος του συστήματος 
u = sin(2*t);
% Οι διαταραχές και ο θόρυβος των μετρήσεων
load w.mat;
load v.mat;
x0 = [0 ; 0];
[out,~,statesB] = lsim(sysmd,[u  w  v],t,x0);
out = out' ;
% Διάνυσμα μετρήσεων
z = [ out(1,:) ; out(2,:)];
%%
% Ορισμός αρχικού πίνακα συνδιασπορας σφαλματος 
% Επιλογή ορθοκανονικού πίνακα ώστε να υπάρχουν
% όλες οι ορίζουσες Di
L = orth(rand(5));
% Ο πίνακας αυτός εγγυάται ότι θα είναι θετικά
% ορισμένος
M = diag(abs(rand(5, 1)) + 0.5);
% Δημιουργία Συμμετρικού πίνακα
P0 = L*M*L';
% Έλεγχος ιδιοτιμών
eig_P = eig(P0);
% Έλεγχος οριζουσών Di
Det_i = [P0(1,1) det(P0(1:2,1:2)) det(P0(1:3,1:3))...
           det(P0(1:4,1:4)) det(P0)];
%%
% Question 3
% Πίνακας συμμεταβλητότητας σφαλμάτων για το
% υποσύστημα που μας ενδιαφέρει
P10 = [P0(4,4) P0(4,5);
       P0(5,4) P0(5,5)];
% Αντίστοιχα για τις καταστάσεις που θεωρούμε 
% αδιάφορες
P20 = [ P0(1,1:3); P0(2,1:3); P0(3,1:3)];
%
S0 = [ P0(4,1:3); P0(5,1:3)];
% Προσδιορισμός παραμέτρων για την υλοποίηση 
% του φίλτρου
A2 = [sysBd.A(1,1:3); sysBd.A(2,1:3); sysBd.A(1,1:3)]; 
H2 = [sysBd.C(1,1:3) ; sysBd.C(2,1:3)];
q2 = 2.5;
Q2 = q2*eye(3);
%%
% Δέσμευση αποθηκευτικού χώρου για μεγέθη που
% θέλουμε να υπολογίσουμε
innovationB = zeros(2,length(t));
x_eB = zeros(2,length(t));
% Αρχικές Εκτιμήσεις για την κατάσταση
x = x0;
% Αρχικές τιμές πινάκων συμμεταβλητότητας σφάλματος
%P1 = P10; P2 = P20; S = S0;
P1 = P10; P2 = P20; S = S0;
% Υλοποίηση Αλγορίθμου του φίλτρου Schmidt-Kalman
for i=1:length(t)
    % Πρόβλεψη
    % x[n](a priori) = A1*x[n-1]+B1*u[n]
    x = A1*x+B1*u(i,:);
    % P1[n](a priori) = A1*P1[n-1]*A1'+Q1
    P1 = A1*P1*A1'+Q1;
    % S[n](a priori) = A1*S[n-1]*A2'
    S = A1*S*A2';
    % P2[n](a priori) = A2*P2[n-1]*A2'+Q2
    P2 = A2*P2*A2'+Q2;
    
    % Διόρθωση
    % an = H1*P1[n](a priori)*H1' + H1*S[n](a priori)*H2'
    %      + H2*S'[n](a priori)*H1' + H2*P2[n](a priori)*H2' + R
    an = H1*P1*H1' + H1*S*H2'+ H2*S'*H1' + H2*P2*H2' + R;
    % K[n] = (P1[n](a priori)*H1' + S[n](a priori)*H2')*inv(an)
    K = (P1*H1'+S*H2')/an;
    % Υπολογισμός καινοτομιών
    innovationB(:,i) = z(:,i) - H1*x ;
    % x[n](a posteriori) = x[n](a priori)+K[n]*(z[n]-H1*x[n](a priori)
    x = x + K*innovationB(:,i);
    % Αποθήκευση εκτίμησης 
    x_eB(:,i) = x; 
    % P1[n](a posteriori) = (I-K[n]*H1)*P1[n](a priori)
    %                     -K[n]*H2*S'[n](a priori)
    P1 = (eye(2)-K*H1)*P1-K*H2*S';
    % S[n](a posteriori) = (I-K[n]*H1)*S[n](a priori)
    %                      -K[n]*H2*P2[n](a priori)
    S = (eye(2)-K*H1)*S - K*H2*P2;
    % P2[n](a posteriori) = P2[n](a priori)
    P2= P2;
end
%%
% Δέσμευση χώρου για μεγαλύτερη ταχύτητα 
% εκτέλεσης
errorB = zeros(2,length(t));
% Απόδοση της εκτίμησης του φίλτρου Kalman
for j = 1:2
    % Εύρεση σφάλματος μεταξύ των καταστάσεων 
    % και της εκτίμησης αυτών
    errorB(j,:) = statesB(:,j) - (x_eB(j,:))';
end
%%
% Γραφικές Παρστάσεις σφαλμάτων
plot(t,errorB(1,:));
grid minor;
title('Estimation Error of x4');
xlabel('Time(s)');
pause();
close;

plot(t,errorB(2,:));
grid minor;
title('Estimation Error of x5');
xlabel('Time(s)');
pause();
close;
%%
% Μέσος όρος απόδοσης φίλτρου Kalman
% για κάθε κατάσταση
N=length(t);
pB = zeros(2,1);
for j = 1:2
    % Κανονικοποίηση ως προς την μέγιστη απόκλιση
    errorB(j,:) = errorB(j,:)/max(abs(errorB(j,:)));
    % Εύρεση νόρμας των διανυσμάτων σφάλματος
    pB(j) = norm(errorB(j,:));
end
% Τύπος RMSE 
performanceB = 1-sqrt(sum(pB)/N)
%%
% Question 4
load x_eA.mat
plot(t,x_eA(4,:),'-r');
hold on;
plot(t,x_eB(1,:),'-b');
grid minor;
title('Estimations of x4');
legend('Kalman Filter','Reduced Kalman Filter');
xlabel('Time(s)');
pause();
close;

plot(t,x_eA(5,:),'-r');
hold on;
plot(t,x_eB(2,:),'-b');
grid minor;
title('Estimations of x5');
legend('Kalman Filter','Reduced Kalman Filter');
xlabel('Time(s)');
pause();
close;
%%
% Απεικονίσεις Καινοτομιών
% Κανονικοποίηση των καινοτομιών ώστε 
% να έχουν εύρος απο [-1 ,1] και να φανεί καλύτερα 
% η μέση τιμή 
innovationB1_n = innovationB(1,:)/max(abs(innovationB(1,:)));
innovationB2_n = innovationB(2,:)/max(abs(innovationB(2,:)));
% Μέσες τιμές καινοτομών 
means_in = [ sum(innovationB1_n)/N ; 
             sum(innovationB2_n)/N ];
% Μέσες τιμές λευκών θορύβων για συγκριση
mean_w = sum(w)/length(w);
mean_v = sum(v)/length(v);
mean_noise = [ mean_w mean_v];
%%
% Απεικόνιση ιστογραμμάτων για να εξετάσουμε αν 
% ικανοποιούν Gaussian κατανομή
histogram(innovationB1_n)
title('Innovation of x4');
pause();
close;
histogram(innovationB2_n)
title('Innovation of x5');
pause();
close;


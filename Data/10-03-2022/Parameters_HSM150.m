clear;

%% Parametre HSM150
Tvz = 1e-3;
Tgm = 1e-3;

% J0 = 1.2e-4;
% B0 = 7.03e-5;

J0 = 1.2153e-4;
B0 = 6.9498e-5;

Mn = 0.39;
Mmax = 0.75*Mn;

N = 2500;

KFI = 0.75;
Cu = 0.0458;
Tn = 0.001;

% Doba simulacie 
Tsim = 20;
% Tsim = 30;

%% Parametre simulovanej sustavy
% Obmedzenia
A = [100; 150];
w_ref = [1; 1];

Dq_max = A.*w_ref;
D2q_max = A.*w_ref.^2;
v_lim = [200; 200];

% Norma vlnoveho cisla
kN = 2*pi/100;

% Matica J
% JN11 = 1.75;
% JN12 = 0.25;
% J_max(1, 1) = J0*(JN11 + JN12);
% J_max(1, 2) = 0;
% J_max(2, 1) = J_max(1, 2);
% J_max(2, 2) = J0;
JN11 = 2;
JN12 = 1;
JN13 = 2;
JN21 = 1;
JN22 = 1;
J_max(1, 1) = J0*(JN11 + JN12 + JN13);
J_max(1, 2) = J0*(JN21 + JN22);
J_max(2, 1) = J_max(1, 2);
J_max(2, 2) = J0;

% Coriolisove sily
CN1 = 0.07/(2*v_lim(1)*v_lim(2) + v_lim(2)^2);
CN2 = 0.07/v_lim(1)^2;
C_max(1, 1) = CN1*(2*v_lim(1)*v_lim(2) + v_lim(2)^2);
C_max(2, 1) = CN2*v_lim(1)^2;

% Gravitacne sily
GN11 = 0.05;
GN12 = 0.03;
GN2 = 0.03;
G_max(1, 1) = GN11 + GN12;
G_max(2, 1) = GN2;

% Nahradne parametre systemu
B = [B0 0; 0 B0];
K = 1./diag(B);
T = [J_max(1, 1); J_max(2, 2)]./diag(B);

% Obmedzenia momentov sil
Mz_max = [J_max(1, 2); J_max(2, 1)].*D2q_max(2:-1:1) + C_max + G_max;
M_max = [0.75*Mn; 0.75*Mn];

M(1, 1) = M_max(1) + Mz_max(1) + v_lim(1)/K(1);
M(1, 2) = M_max(1) - Mz_max(1) - v_lim(1)/K(1);
M(2, 1) = M_max(2) + Mz_max(2) + v_lim(2)/K(2);
M(2, 2) = M_max(2) - Mz_max(2) - v_lim(2)/K(2);

mu = [
    M(1, 2)/M_max;
    M(2, 2)/M_max
];

% Parametre zakona riadenia
d = [1; 1];

w_max = [
    (-1 + sqrt(1 + 4*K(1)*M(1, 2)*T(1)/A(1)))/(2*T(1));
    (-1 + sqrt(1 + 4*K(2)*M(2, 2)*T(2)/A(2)))/(2*T(2))
];
a = [
    1/(T(1)*(1 - K(1)*M(1, 2)/v_lim(1)*log(1 + v_lim(1)/(K(1)*M(1, 2)))));
    1/(T(2)*(1 - K(2)*M(2, 2)/v_lim(2)*log(1 + v_lim(2)/(K(2)*M(2, 2)))))
];
a = [
    1/(T(1)*(1 - K(1)*M(1, 2)/v_lim(1)/2*log(1 + 2*v_lim(1)/(K(1)*M(1, 2)))));
    1/(T(2)*(1 - K(2)*M(2, 2)/v_lim(2)/2*log(1 + 2*v_lim(2)/(K(2)*M(2, 2)))))
];
a = [
    1/(T(1)*(1 - M(1, 2)/M(1, 1)*log(1 + M(1, 1)/M(1, 2))));
    1/(T(2)*(1 - M(2, 2)/M(2, 1)*log(1 + M(2, 1)/M(2, 2))))
];
% w_max = a;
k = [
    w_max(1)*100/(mu(1)*d(1))*sqrt((w_max(1)^2 + 1/T(1)^2)/(w_max(1)^2 + a(1)^2));
    w_max(2)*100/(mu(2)*d(2))*sqrt((w_max(2)^2 + 1/T(2)^2)/(w_max(2)^2 + a(2)^2))
];

Kp = T.*a.*k./K;
Kd = (T.*(a + k) - 1)./K;

%% Vykreslenie
t = data(:, 1);
q1_ref = data(:, 2);
q2_ref = data(:, 3);
q1 = data(:, 4);
q2 = data(:, 5);
Dq1_ref = data(:, 6);
Dq2_ref = data(:, 7);
Dq1 = data(:, 8);
Dq2 = data(:, 9);
e1 = data(:, 10);
e2 = data(:, 11);
De1 = data(:, 12);
De2 = data(:, 13);
u1 = data(:, 14);
u2 = data(:, 15);
Mm = data(:, 16);
Mz = data(:, 17);
M01 = data(:, 18);
M02 = data(:, 19);
I = data(:, 20);

figure;
    subplot(211);
        hold on;
            plot(t, q1, 'r');
            plot(t, q2, 'b');
            plot(t, q1_ref, 'k:');
            plot(t, q2_ref, 'k:');
        hold off;
        title('Poloha');

    subplot(212);
        hold on;
            plot(t, Dq1, 'r');
            plot(t, Dq2, 'b');
            plot(t, Dq1_ref, 'k:');
            plot(t, Dq2_ref, 'k:');
        hold off;
        title('Rychlost');

figure;
    hold on;
        plot(t, e1*100/A(1), 'r');
        plot(t, e2*100/A(2), 'b');
    hold off;
    title('Odchylka');

figure;
    subplot(311);
        hold on;
            plot(t, u1, 'r');
            plot(t, u2, 'b');
        hold off;
        title('Akcny zasah');
    subplot(312);
        plot(t, Mm, 'b');
        title('Akcny zasah - meranie');
    subplot(313);
        plot(t, Mz, 'b');
        title('Porucha');
% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Mouton Laudin, Alfonso
% Tp N� 1 - Caso de estudio 1 - 
%
%   Inciso 5  
%   A partir de las curvas de mediciones de las variables graficadas en la Fig. 1-3, se requiere 
%   obtener el modelo del sistema considerando como entrada un escal�n de 12V, como salida a la 
%   velocidad angular, y al torque de carga TL aplicado una perturbaci�n. En el archivo 
%   Curvas_Medidas_Motor.xls est�n las mediciones, en la primer hoja los valores y en la segunda 
%   los nombres. Se requiere obtener el modelo din�mico, para establecer las constantes del modelo 
%   (1-5) (1-6).  
%%

clc; clear all; close all;

%Lectura de los datos del xls
Datos = xlsread('Curvas_Medidas_Motor_2024.xls',1);
t = Datos(1:end,1);
Wm = Datos(1:end,2);
I = Datos(1:end,3);
Va = Datos(1:end,4);
Ei = 12;
deltaT = Datos(1,1);
delayT = Datos(702,1);

%Visualizaci�n de los Datos del Excel
figure(1);
subplot(3,1,1)
plot(t,Wm); grid;
ylabel('[rad/s]'); xlabel('[s]'); title('Velocidad Angular'); ylim([0 250]);
subplot(3,1,2)
plot(t,I); grid;
ylabel('[A]'); xlabel('[s]'); title('Corriente');
subplot(3,1,3)
plot(t,Va); grid;
ylabel('[V]'); xlabel('[s]'); title('Tensi�n Armadura');

%Acercando en el transitorio
figure(2);
zoom = [695 715];
subplot(3,1,1)
plot(t(zoom(1):zoom(2)) ,Wm(zoom(1):zoom(2))); grid;
ylabel('[rad/s]'); xlabel('[s]'); title('Velocidad Angular'); ylim([0 250]);
subplot(3,1,2)
plot(t(zoom(1):zoom(2)) ,I(zoom(1):zoom(2))); grid;
ylabel('[A]'); xlabel('[s]'); title('Corriente');
subplot(3,1,3)
plot(t(zoom(1):zoom(2)) ,Va(zoom(1):zoom(2))); grid;
ylabel('[V]'); xlabel('[s]'); title('Tensi�n Armadura');

%M�todo de chen
%Elecci�n de 3 puntos equisdistantes:
inicio_W = 703;
delta_W = 2;
t1_W = Datos(inicio_W,1); W1 = Datos(inicio_W,2);
t2_W = Datos(inicio_W + delta_W,1); W2 = Datos(inicio_W + delta_W,2);
t3_W = Datos(inicio_W + 2*delta_W,1); W3 = Datos(inicio_W + 2*delta_W,2);

inicio_I = 703;
delta_I = 1;
t1_I = Datos(inicio_I,1); I1 = Datos(inicio_I,3);
t2_I = Datos(inicio_I + delta_I,1); I2 = Datos(inicio_I + delta_I,3);
t3_I = Datos(inicio_I + 2*delta_I,1); I3 = Datos(inicio_I + 2*delta_I,3);

%Visualizaci�n de los puntos en el 'Zoom'
figure(3);
subplot(2,1,1)
plot(t(zoom(1):zoom(2)) ,Wm(zoom(1):zoom(2)), [t1_W t2_W t3_W],[W1 ,W2 ,W3],'+'); grid;
ylabel('[rad/s]'); xlabel('[s]'); title('Velocidad Angular'); ylim([0 250]);
subplot(2,1,2)
plot( t(zoom(1):zoom(2)) ,I(zoom(1):zoom(2)) , [t1_I t2_I t3_I],[I1 ,I2 ,I3],'+'); grid;
ylabel('[A]'); xlabel('[s]'); title('Corriente');

%Obtengo las ganancias: (Y las normalizo)
k_W = Datos(710,2)/Ei; k_I = Datos(710,3)/Ei;

k1_W = (W1/(Ei*k_W))-1; k2_W = (W2/(Ei*k_W))-1; k3_W = (W3/(Ei*k_W))-1;
k1_I = (I1/(Ei*k_I))-1; k2_I = (I2/(Ei*k_I))-1; k3_I = (I3/(Ei*k_I))-1;

%Aplico las ec 23

%OMEGA_MOT
b_W = 4*k1_W^3*k3_W - 3*k1_W^2*k2_W^2 - 4*k2_W^2 + k3_W^2 + 6*k1_W*k2_W*k3_W;
alpha1_W = (k1_W*k2_W + k3_W - sqrt(b_W))/(2*(k1_W^2 + k2_W));
alpha2_W = (k1_W*k2_W + k3_W + sqrt(b_W))/(2*(k1_W^2 + k2_W));
%beta_W = (2*k1_W^3 + 3*k1_W*k2_W + k3_W - sqrt(b_W))/(sqrt(b_W));
beta_W = (k1_W + alpha2_W)/(alpha1_W - alpha2_W);

T1_W = -(t1_W - delayT)/log(alpha1_W);
T2_W = -(t1_W - delayT)/log(alpha2_W);
T3_W = beta_W*(T1_W-T2_W) + T1_W;

G_W = tf(k_W*[T3_W 1],conv([T1_W 1],[T2_W 1]))

%CORRIENTE
b_I = 4*k1_I^3*k3_I - 3*k1_I^2*k2_I^2 - 4*k2_I^2 + k3_I^2 + 6*k1_I*k2_I*k3_I;
alpha1_I = (k1_I*k2_I + k3_I - sqrt(b_I))/(2*(k1_I^2 + k2_I));
alpha2_I = (k1_I*k2_I + k3_I + sqrt(b_I))/(2*(k1_I^2 + k2_I));
%beta_I = (2*k1_I^3 + 3*k1_I*k2_I + k3_I - sqrt(b_I))/(sqrt(b_I));
beta_I = (k1_I + alpha2_I)/(alpha1_I - alpha2_I);

T1_I = -(t1_I - delayT)/log(alpha1_I);
T2_I = -(t1_I - delayT)/log(alpha2_I);
T3_I = beta_I*(T1_I - T2_I) + T1_I;

G_I = tf(k_I*[T3_I 1],conv([T1_I 1],[T2_I 1]))

%Se compara la FT obtenida con los datos del excel:
[Wobt, tobt] = lsim(G_W, Va(1:1500), t(1:1500));
[Iobt, tobt] = lsim(G_I, Va(1:1500), t(1:1500));

figure(4)
subplot(4,1,1)
plot(tobt,Wm(1:1500));
title('Velocidad Angular Excel'); xlabel('Tiempo [s]'); ylabel('[Rad/s]'); ylim([0 250]);
subplot(4,1,2)
plot(tobt,Wobt);
title('Velocidad Angular Obtenida'); xlabel('Tiempo [s]'); ylabel('[Rad/s]'); ylim([0 250]);
subplot(4,1,3)
plot(tobt,I(1:1500));
title('Corriente Excel'); xlabel('Tiempo [s]'); ylabel('[I]');
subplot(4,1,4)
plot(tobt,Iobt);
title('Corriente Obtenida'); xlabel('Tiempo [s]'); ylabel('[I]');

%con la FT GENERAL obtenemos los PAR�METROS del motor:
ki = 1/(G_W.num{1}(3))
km = ki
R = Ei/max(Iobt)
J = G_I.num{1}(2)
B = G_I.num{1}(3)
L = G_I.den{1}(1)*ki*km/J


clear all; clc; close all;

%Sistema de dos variables de estado
R = 47; L = 1e-6; C=100e-9; %Valores asignados Para el RLC
Ei = 12;                    %Tensi�n de entrada

A = [0 1/C ; -1/L -R/L]     %Matriz de Estados
B = [0 1/L]'                %Matriz de entrada a Estado
C1 = [0 R]                  %Matriz de estado a Salida (VR)
C2 = [1 0]                  %Matriz de estado a Salida (VC)
C3 = [0 1]                  %Matriz de estado a Salida (I)
D = 0                       %Matriz de transmisi�n directa

x0 = [0 0]';                %Condici�nes iniciales
Sys1 = ss(A,B,C1,D)         %Sistema (Salida VR)
Sys2 = ss(A,B,C2,D)         %Sistema (Salida VC)
Sys3 = ss(A,B,C3,D)         %Sistema (Salida I)

%Simulaci�n Del Sistema
Ti = 2e-3;                       %Per�odo de la se�al de entrada
N = 1;                           %N�mero de per�odos a graficar
t = 0:Ti/3000:Ti*N;              %Vector tiempo
u = Ei*square(2*pi/Ti*(t))       %Se�al Cuadrada con per�odo 2mS


subplot(4,1,1)
lsim(Sys1, u, t)
title('Resistencia')
xlabel('Tiempo [s]');
ylabel('Tensi�n [V]');
ylim([-30 30])

subplot(4,1,2)
lsim(Sys2, u, t)
title('Capacitor')
xlabel('Tiempo');
ylim([-20 20])
ylabel('Tensi�n [V]');

subplot(4,1,3)
lsim(Sys3, u, t)
title('Corriente')
xlabel('Tiempo');
ylim([-0.8 0.8])
ylabel('Corriente [A]');

subplot(4,1,4)
plot(t, u)
title('Entrada u(t)')
xlabel('Tiempo');
ylim([-13 13])
ylabel('Tensi�n [V]');
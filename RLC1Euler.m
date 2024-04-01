% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Mouton Laudin, Alfonso
% Tp N� 1 - Caso de estudio 1 - 
%   Inciso 1  
%   Asignar valores a R=47ohm, L=1uHy, y C=100nF. Obtener simulaciones que permitan 
%   estudiar la din�mica del sistema, con una entrada de tensi�n escal�n de 12V, que cada 
%   1ms cambia de signo.
%   
%%
clear all; clc; close all;

%Sistema de dos variables de estado
R = 47; L = 1e-6; C = 100e-9;           %Valores asignados Para el RLC
Ei = 12;                                %Tensi�n de entrada

A = [0 1/C ; -1/L -R/L]                 %Matriz de Estados
B = [0 1/L]'                            %Matriz de entrada a Estado
C = [0 R]                               %Matriz de estado a Salida (VR)
D = 0                                   %Matriz de transmisi�n directa

[num, den] = ss2tf(A,B,C,D);            %Sistema (Salida VR)
G = tf(num, den)
polos = pole(G)

%--------------------------------------------------------------------------
tL = (log(0.05)/(polos(2)))*3            %Respuesta m�s lenta del sistema
tR = (log(0.95)/(polos(1)))/3            %Respuesta m�s r�pida del sistema

%CONDICIONES INICIALES:
Y(1) = 0;
I(1) = 0;
VC(1) = 0;
x = [VC(1) I(1)]';
xop = [0 0]';

%Simulaci�n mediante Euler
tF = 2e-3;                               %Tiempo final de la simulaci�n
N = floor(tF/tR)                         %N�mero de puntos a simular
t = 0:tR:N*tR;                           %Vector tiempo
u = Ei*square(2*pi/2e-3*(t));            %Se�al Cuadrada con per�odo 2mS

for i=1: N
    xp = A*(x-xop) + B*u(i);
    x = x + xp*tR;
    aux = C*x;
    Y(i+1) = aux(1);
    I(i+1) = x(2);
    VC(i+1) = x(1);
end

%GR�FICAS
subplot(4,1,1)
plot(t, Y);
title('Resistencia')
xlabel('Tiempo [s]');
ylabel('Tensi�n [V]');
ylim([-30 30])

subplot(4,1,2)
plot(t, VC);
title('Capacitor')
xlabel('Tiempo');
ylim([-20 20])
ylabel('Tensi�n [V]');

subplot(4,1,3)
plot(t, I);
title('Corriente')
xlabel('Tiempo');
ylim([-0.8 0.8])
ylabel('Corriente [A]');

subplot(4,1,4)
plot(t, u);
title('Entrada u(t)')
xlabel('Tiempo');
ylim([-13 13])
ylabel('Tensi�n [V]');
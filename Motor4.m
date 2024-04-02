% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Mouton Laudin, Alfonso
% Tp N° 1 - Caso de estudio 1 - 
%
%   Inciso 4  
%   Obtener el torque máximo que puede soportar el motor modelado mediante las Ecs. (1-5) 
%   (1-6) y (1-7) cuando se lo alimenta con 12V, graficando para 5 segundos de tiempo la velocidad 
% 	angular y corriente ia para establecer su valor máximo como para dimensionar dispositivos 
%   electrónicos. 
%%

clc; clear all; close all;

Laa = 366e-6; J = 5e-9; Ra = 55.6; Bm = 0; Ki = 6.49e-3; Km = 6.53e-3;     %Parámetros del motor.
Ei = 12;                                                                   %Tensión de entrada

%Variables de estados, x1 = ia; x2 = wr; x3 = theta; u = Va, Tl
%Sistema M I S O

A = [-Ra/Laa -Km/Laa 0 ; Ki/J -Bm/J 0 ; 0 1 0]
B = [1/Laa 0 ; 0 -1/J ; 0 0]
C = [0 1 0]                                                                %Salida omega_R
D = [0 0]

%Datos para la simulación
deltaT = 10e-6;                                                           %Paso para Euler
tF = 5;                                                                    %Tiempo final para la simulación
N = floor(tF/deltaT)                                                       %N puntos a simular

%Condiciones Iniciles
ia(1) = 0;
wr(1) = 0;
theta(1) = 0;
x = [ia(1) wr(1) theta(1)]';
xop = [0 0 0]';

%Simulación
t = 0:deltaT:N*deltaT;
u = Ei*ones(1, length(t));                                                 %Entrada de tensión Va
TORQUE_MAX = 0.0014
%Perturbación del tipo escalon
    tL = TORQUE_MAX*heaviside(t - 2);
%Perturbación del tipo rampa
    %tL = TORQUE_MAX*(t-2);                                                     
    %tL(tL>TORQUE_MAX) = TORQUE_MAX;
    %tL(tL<0) = 0;
    %plot(t, tL); grid;

for i=1 : N
    xp = A*(x-xop) + B*[u(i) tL(i)]';
    x = x + xp*deltaT;                                                     %Euler                                                         
    wr(i+1) = C*x;
    ia(i+1) = x(1);
    theta(i+1) = x(3);
end

%GRÁFICAS
subplot(4,1,1)
plot(t, wr); grid;
title('Motor')
xlabel('Tiempo [s]');
ylabel('Velocidad angular [Rad/s]');

subplot(4,1,2)
plot(t, ia); grid;
title('Armadura')
xlabel('Tiempo [s]');
ylabel('Corriente [A]');

subplot(4,1,3)
plot(t, theta);
title('Motor')
xlabel('Tiempo [s]');
ylabel('Angulo [Rad]');

subplot(4,1,4)
plot(t, tL*100*100, t, u);
title('Entradas')
xlabel('Tiempo [s]');
ylabel('[m.N.cm][Volts]');
legend('Perturbación','Tensión Armadura')
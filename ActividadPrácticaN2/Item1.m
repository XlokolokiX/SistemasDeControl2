%% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Mouton Laudin, Alfonso
% Tp N° 2 - ITEM 1 - 
%
%   Implementar un sistema en variables de estado que controle el ángulo del motor, para 
%   consignas de pi/2 y –pi/2 cambiando cada 5 segundos y que el TL es el descripto en la planilla de datos 
%   comparando el desempeño con el obtenido con el PID digital del TP Nº1. Hallar el valor de 
%   integración Euler adecuado.
%
%   OBJETIVO: acelerar la dinámica del controlador verificando el resultado con las curvas del archivo xlsx 
%   adjunto.
% 
%   -Evitando que la tensión supere los 24Volts en valor absoluto, especificar el tiempo de muestreo 
%    necesario para el controlador cumpla el objetivo. 
%   -Asumiendo que no puede medirse directamente la corriente, pero sí la velocidad y el ángulo, 
%    proponer un controlador que logre el objetivo. 
%   -Determinar el efecto de la nolinealidad en la acción de control, descripta en la Fig. 2, y verificar cuál 
%    es el máximo valor admisible de ésa no linealidad. 
%%
clc; clear all; close all;

%MEDICIONES EXCEL:
Datos = xlsread('Curvas_Medidas_Motor_2024.xls',1);
t_E = Datos(1:end,1);
Wr_E = Datos(1:end,2);
I_E = Datos(1:end,3);
Va_E = Datos(1:end,4);
Tl_E = Datos(1:end,5);
Ei_E = 12;

deltaT = Datos(1,1);
delayT = Datos(702,1);

max(Wr_E)*2
max(I_E)*2

%%
%%MODELADO DEL MOTOR LAZO ABIERTO (Contínuo)
%Parámetros del motor obtenidos anteriormente por CHEN
Ki = 0.01162; J = 2.0628e-9; Bm = 0; Laa = 7.0274e-4; Ra = 28.13; Km = 0.0605;

A_c = [-Ra/Laa -Km/Laa 0;     %x1: Corriente
        Ki/J  -Bm/J    0;     %x2: Velocidad Angular
        0       1      0]     %x3: Angulo

B_c = [1/Laa 0;
        0 -1/J;
        0    0]

C_c = [0  0  1]               %Ángulo

D_c = [0 0]

G = ss(A_c, B_c, C_c, D_c)

%Determinación del los tiempos de la dinámica del sistema contínuo
root_c = eig(A_c);
tR = log(0.95)/ real(root_c(2))
tL = log(0.05)/ real(root_c(2))

%%
%SISTEMA DISCRETIZADO
Tm = (tR/4)                 %Tiempo de muestreo

G_d = c2d(G ,Tm ,'zoh');    %Discretización mediante un Retentor de Orden Cero
A_d = G_d.a
B_d = G_d.b;
B_d = B_d(:,1)              %Ya que sólo me interesa controlar el ángulo
C_d = G_d.c
D_d = G_d.d
%%
%Prueba de Controlabilidad - Alcanzabilidad
Ma = [B_d  A_d*B_d  A_d^2*B_d];
rank(Ma)
Mc = [B_d  A_d*B_d  A_d^2*B_d A_d^3];
rank(Mc)

%%
%SISTEMA AMPLIADO
AA = [A_d , zeros(3,1) ; -C_d(1,:)*A_d, 1]
BA = [B_d ; -C_d(1,:)*B_d]

%LQR CONTROLADOR 
%Se dimensiona a prueba y error
%Ponderaciones mayores -> Se penalizará más el mismo (Indica un elevado costo)

%d_c = [1.3744 6.25e-6 0.40528 0.1];        %Corriente, velocidad , Ángulo, psita
d_c = [1.3737 6.3477e-6 0.40528 1];        %Corriente, velocidad , Ángulo, psita
Q_c = diag(d_c);                           %Q actua sobre las VARIABLES DE ESTADO
R_c = 3e4;                                 %R actua sobre la ACCION DE CONTROL

KLQR_C = dlqr(AA, BA, Q_c, R_c);
K = KLQR_C(1:3);
Ki = -KLQR_C(4);

%%
%Sistema DUAL
Ao = A_d'
Bo = C_d'
Co = B_d'

%LQR OBSERVADOR
d_o = [10 10 .09e0];                        %Corriente, velocidad, Ángulo
Q_o = diag(d_o);                            %Matriz de covarianza del ruido del proceso (incertidumbre de variables -> Valores elevados)
R_o = 1e1                                   %Matriz de covarianza del ruido de la medición (incertidumbre en las mediciones -> Valores elevados)
K_o = (dlqr(Ao,Bo,Q_o,R_o))'

%%
%SIMULACIÓN
Bloques_Sist = imread('Bloques_Motor.png');
figure(1);
imshow(Bloques_Sist);

Tf = 2;                                 %Tiempo final para la simulación
Ti = Tm                                 %Tiempo de Integración
N = floor(Tf/Ti)                        %Numero de pasos para la simulación
death_zone = 1;                         %Zona muerta del Actuador

t = 0:Ti:N*Ti;
Wr = zeros(1, N+1);
Theta = zeros(1, N+1);
Ia = zeros(1, N+1);
u = zeros(1, N+1);

%Referencia que cambia cada 5seg
%Ref = (pi/2)*square(t*2*pi/5); 
Ref = (pi/2)*square(t*2*pi/0.5); 
%Vector perturbación
Tl = max(Tl_E)*square(t*2*pi/0.3);
Tl(Tl<0)=0;

%Condiciones_Iniciales
Wr(1) = 0;
Theta(1) = 0;
Ia(1) = 0;
x = [Ia(1) Wr(1) Theta(1)]';
x_obs = [0 0 0]';
psita = 0;

for i=1: N
    
    y = C_d*x;
    y_obs = C_d*x_obs;
    
    psita_p = Ref(i) - C_d*x;       %Err Integrador
    psita = psita + psita_p;
    
    %Control
    %u(i) = -K*x + Ki*psita;          %Sin Observador
    u(i) = -K*x_obs + Ki*psita;      %Con Observador
    
    %Comportamiento No Lineal del Actuador
    if(abs(u) < death_zone)
        u(i) = 0;
    else
        u(i) = sign(u(i))*( abs(u(i))-death_zone );
    end
    
    %SISTEMA
    Ia_p = -(Ra/Laa)*Ia(i) -(Km/Laa)*Wr(i) + (1/Laa)*u(i);
    Ia(i+1) = Ia(i) + Ia_p*Ti;
    Wr_p = (Ki/J)*Ia(i) - (Bm/J)*Wr(i) - (1/J)*Tl(i);
    Wr(i+1) = Wr(i) + Ti*Wr_p;
    Theta(i+1) = Theta(i) + Ti*Wr(i);
    
    x_obs = A_d*x_obs + B_d*u(i) + K_o*(y - y_obs);
    x = [Ia(i+1) Wr(i+1) Theta(i+1)]';
end

%%
%GRÁFICAS

figure(2);
subplot(3,2,1);hold on;
plot(t ,Theta,t,Ref); title('Ángulo  \phi [rad]'); grid on; hold on;
subplot(3,2,2);hold on;
plot(t ,Wr); title('Velocidad Angular  \omega [rad/s]');grid on;hold on; 
subplot(3,2,3);hold on;
plot(t ,Ia); title('Corriente  Ia [A]');grid on;hold on; 
subplot(3,2,4);hold on;
plot(t ,u); title('Acción de control  V [V]');grid on;hold on;
subplot(3,1,3);hold on;
plot(t ,Tl);title('Torque  TL [N/m]');grid on;hold on;

figure(3)
plot(Theta ,Wr); title('Plano de Fases'); grid on;hold on;


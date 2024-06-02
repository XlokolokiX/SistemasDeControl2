%% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Mouton Laudin, Alfonso
% Tp N° 2 - ITEM 2 - 
%
% Para el caso del avión, emplear un tiempo de integración por Euler adecuado y un tiempo de 
%simulación de 70seg. Los parámetros son a=0.07; ?=9; b=5; c=150, hallar un controlador para que 
%los polos de lazo cerrado se ubican en ui=-15+-15j; -0.5+-0.5j, para referencias de 100 y -100 metros
%en altura, ambas con alturas iniciales de -500 y 500. 
%  -Proponer un controlador en tiempo discreto en variables de estado para que el proceso evolucione 
%en los rangos de validez del modelo, es decir donde los ángulos y el valor de la acción de control en 
%valor absoluto son menores a la unidad. 
%  -Asumiendo que no puede medirse el ángulo ?, pero sí el ángulo ? y la altura, proponer un esquema 
%que permita lograr el objetivo de control. 
%  -Establecer el valor del tiempo de muestreo más adecuado para implementar el diseño en un sistema 
%micro controlado. 
%  -Determinar el efecto de la nolinealidad en la acción de control, descripta en la Fig. 4, y verificar cuál 
%es el máximo valor admisible de la nolinealidad. 
%%
clc; clear all; close all;

%MODELADO DEL SISTEMA EN EE. (Contínuo)
a = 0.07; w = 9; b = 5; c = 150;

A_c = [a   -a   0  0;       %x1 = angulo de cabeceo
       -w^2 w^2 0  0;       %x2 = angulo
       0    c   0  0;       %x3 = altura
       0    0   0  1]       %x4 = theta_p
   
B_c =  [0;
        0;
        c;
        0]

C_c = [0   0   1   0]

D_c = [0]

G = ss(A_c, B_c ,C_c , D_c)
%Determinación de los tiempos de dinámica
roots_c = eig(A_c);
tR = abs(log(0.95)/real(roots_c(4)))

%%
%SISTEMA DISCRETIZADO
Tm = tR/4;

G_d = c2d(G, Tm , 'zoh');
A_d = G_d.a
B_d = G_d.d
C_d = G_d.c
D_d = G_d.d

%%
%Prueba de Controlabilidad - Alcanzabilidad
Ma = [B_d  A_d*B_d  A_d^2*B_d A_d^3*B_d];
rank(Ma)
Mc = [B_d  A_d*B_d  A_d^2*B_d A_d^3*B_d A_d^3];
rank(Mc)

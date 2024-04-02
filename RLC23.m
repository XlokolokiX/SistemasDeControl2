% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Mouton Laudin, Alfonso
% Tp N° 1 - Caso de estudio 1 - 
%
%   Inciso 2-3  
%   En el archivo Curvas_Medidas_RLC.xls (datos en la hoja 1 y etiquetas en la hoja 2) 
%   están las series de datos que sirven para deducir los valores de R, L y C del circuito. 
%   Emplear el método de la respuesta al escalón, tomando como salida la tensión en el capacitor. 
%   Una vez determinados los parámetros R, L y C, emplear la serie de corriente desde 
%   0.05seg en adelante para validar el resultado superponiendo las gráficas.
%%
clear all; clc; close all;

%Lectura de los datos del excel
excel = xlsread('Curvas_Medidas_RLC_2024.xls'); 
t = excel(1:end,1);
I = excel(1:end,2);
Vc = excel(1:end,3);
u = excel(1:end,4);

Ei = 12;

%Gráfica de los datos
subplot(3,1,1)
plot(t,u);
title('Entrada')
xlabel('Tiempo [s]');
ylabel('Tensión [V]');
subplot(3,1,2)
plot(t,Vc);
title('Capacitor')
xlabel('Tiempo [s]');
ylabel('Tensión [V]');
subplot(3,1,3)
plot(t,I);
title('Corriente')
xlabel('Tiempo [s]');
ylabel('Corriente [I]');
%Se aprecia un sistema de 2do Orden sin sobrepasamiento


close all;
%Método de Chen
%Identificación de sistemas de 2do Orden basado en la respuesta al escalón
%Elijo 3 Puntos (equisdistantes) para aplicar el método
delta_V = 2
inicio_V = 104
t1_Vc = excel(inicio_V,1); Vc1 = excel(inicio_V,3)/Ei; 
t2_Vc = excel(inicio_V + delta_V,1); Vc2 = excel(inicio_V + delta_V,3)/Ei;
t3_Vc = excel(inicio_V + delta_V*2,1); Vc3 = excel(inicio_V + delta_V*2,3)/Ei;

delta_I = 1
inicio_I = 102
t1_I = excel(inicio_I,1); I1 = excel(inicio_I,2)/Ei;
t2_I = excel(inicio_I + delta_I,1); I2 = excel(inicio_I + delta_I,2)/Ei;
t3_I = excel(inicio_I + delta_I*2,1); I3 = excel(inicio_I + delta_I*2,2)/Ei;

%Ganancias
K_Vc = excel(500,3)/Ei;     K_I = excel(498,2)/Ei;

%Defino las k correspondientes para los puntos tomados
k1_V = (Vc1/K_Vc) - 1;   k1_I= (I1/K_I) - 1;
k2_V = (Vc2/K_Vc) - 1;   k2_I= (I2/K_I) - 1;
k3_V = (Vc3/K_Vc) - 1;   k3_I= (I3/K_I) - 1;

%Tomando la EC23 del método PARA VC
b_V = 4*k1_V^3*k3_V - 3*k1_V^2*k2_V^2 - 4*k2_V^3 + k3_V^2 + 6*k1_V*k2_V*k3_V
alpha1_V = (k1_V*k2_V + k3_V - sqrt(b_V))/(2*(k1_V^2 + k2_V))
alpha2_V = (k1_V*k2_V + k3_V + sqrt(b_V))/(2*(k1_V^2 + k2_V))
%beta_V = (2*k1_V^3 + 3*k1_V*k2_V + k3_V - sqrt(b_V))/(sqrt(b_V))
beta_V = (k1_V+alpha2_V)/(alpha1_V-alpha2_V)

T1_V = -(t1_Vc - 0.01)/log(alpha1_V)
T2_V = -(t1_Vc - 0.01)/log(alpha2_V)
T3_V = (beta_V*(T1_V - T2_V)) + T1_V

G_V = tf(K_Vc*[T3_V 1], conv([T1_V 1],[T2_V 1]) )

%Tomando la EC23 del método PARA I
b_I = 4*k1_I^3*k3_I - 3*k1_I^2*k2_I^2 - 4*k2_I^3 + k3_I^2 + 6*k1_I*k2_I*k3_I
alpha1_I = (k1_I*k2_I + k3_I - sqrt(b_I))/(2*(k1_I^2 + k2_I))
alpha2_I = (k1_I*k2_I + k3_I + sqrt(b_I))/(2*(k1_I^2 + k2_I))
%beta_I = (2*k1_I^3 + 3*k1_I*k2_I + k3_I - sqrt(b_I))/(sqrt(b_I))
beta_I = (k1_I+alpha2_I)/(alpha1_I-alpha2_I)

T1_I = -(t1_I - 0.01)/log(alpha1_I)
T2_I = -(t1_I - 0.01)/log(alpha2_I)
T3_I = (beta_I*(T1_I - T2_I)) + T1_I

G_I = tf( K_I*[T3_I 1], conv([T1_I 1],[T2_I 1]) )

%Se compara la FT obtenida con la del excel:
[Vobt, tobt] = lsim(G_V, u, t);
[Iobt, tobt] = lsim(G_I, u, t);

subplot(3,1,1)
plot(t,u);
title('Entrada')
xlabel('Tiempo [s]');
ylabel('Tensión [V]');
subplot(3,1,2)
plot(t,Vc,t,Vobt);
title('Capacitor')
xlabel('Tiempo [s]');
ylabel('Tensión [V]');
subplot(3,1,3)
plot(t,I,t,Iobt);
title('Corriente')
xlabel('Tiempo [s]');
ylabel('Corriente [I]');

%De la función de transferencia sacamos
C = G_I.num{1}(2)
L = (G_I.den{1}(1))/C
R = (G_I.den{1}(2))/C

%Matrices
A = [-R/L -1/L; 1/C 0];
B = [1/L; 0];
C = [1 0];
D = 0;

%Definicion de la ecuación de estado y de salida (salida de corriente)
G1 = ss(A,B,C,D);
[yout,yt] = lsim(G1,(u),t);

figure(6)
plot(yt, yout,'b');grid on; hold on;
plot( valores(:,1), valores(:,2),'r' ); title('Comparación de corriente');
legend({'i(t) aproximada con valores RLC calculados','i(t) de excel'},'Location','southeast')

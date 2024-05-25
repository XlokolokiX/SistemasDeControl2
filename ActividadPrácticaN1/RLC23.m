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

%LECTURA EXCEL-------------------------------------------------------------
excel = xlsread('Curvas_Medidas_RLC_2024.xls'); 
t = excel(1:end,1);
I = excel(1:end,2);
Vc = excel(1:end,3);
u = excel(1:end,4);
Ei = max(u);
%--------------------------------------------------------------------------
%GRÁFICA DEL EXCEL---------------------------------------------------------
figure(1)
subplot(3,1,1)
plot(t,u);
title('Entrada'); xlabel('Tiempo [s]'); ylabel('Tensión [V]');
subplot(3,1,2)
plot(t,Vc);
title('Capacitor'); xlabel('Tiempo [s]'); ylabel('Tensión [V]'); 
subplot(3,1,3)
plot(t,I);
title('Corriente'); xlabel('Tiempo [s]'); ylabel('Corriente [I]');
%--------------------------------------------------------------------------
%MÉTODO DE CHEN
%ALGORITMO DE BÚSQUEDA DE LOS PUNTOS---------------------------------------
tI_Inicio = 0.0105;
tI_Delta = 0.0005;
tVc_Inicio = 0.011;
tVc_Delta = 0.001;
[~ , Ind1_I] = min(abs(t-tI_Inicio))
[~ , Ind2_I] = min(abs(t-(tI_Inicio + tI_Delta)))
[~ , Ind3_I] = min(abs(t-(tI_Inicio + 2*tI_Delta)))
[~ , Ind1_Vc] = min(abs(t-tVc_Inicio))
[~ , Ind2_Vc] = min(abs(t-(tVc_Inicio + tVc_Delta)))
[~ , Ind3_Vc] = min(abs(t-(tVc_Inicio + 2*tVc_Delta)))
t1_Vc = excel(Ind1_Vc,1); Vc1 = excel(Ind1_Vc,3);
t2_Vc = excel(Ind2_Vc,1); Vc2 = excel(Ind2_Vc,3);
t3_Vc = excel(Ind3_Vc,1); Vc3 = excel(Ind3_Vc,3);
t1_I = excel(Ind1_I,1); I1 = excel(Ind1_I,2);
t2_I = excel(Ind2_I,1); I2 = excel(Ind2_I,2);
t3_I = excel(Ind3_I,1); I3 = excel(Ind3_I,2);

%Visualización de los puntos
figure(2);
zoom = [90 150];
subplot(2,1,1)
plot(t(zoom(1):zoom(2)) ,I(zoom(1):zoom(2)), [t1_I t2_I t3_I],[I1 ,I2 ,I3],'+'); grid;
title('Corriente'); xlabel('Tiempo [s]'); ylabel('Corriente [I]');
subplot(2,1,2)
plot( t(zoom(1):zoom(2)) ,Vc(zoom(1):zoom(2)) , [t1_Vc t2_Vc t3_Vc],[Vc1 ,Vc2 ,Vc3],'+'); grid;
title('Capacitor'); xlabel('Tiempo [s]'); ylabel('Tensión [V]'); ylim([0 11]);
%--------------------------------------------------------------------------
%GANANCIAS-----------------------------------------------------------------
K_Vc = excel(500,3)/Ei;     K_I = excel(498,2)/Ei;
k1_V = (Vc1/(Ei*K_Vc)) - 1;   k1_I= (I1/(Ei*K_I)) - 1;
k2_V = (Vc2/(Ei*K_Vc)) - 1;   k2_I= (I2/(Ei*K_I)) - 1;
k3_V = (Vc3/(Ei*K_Vc)) - 1;   k3_I= (I3/(Ei*K_I)) - 1;
%--------------------------------------------------------------------------
%RESOLVIENDO (EC23) PARA VC
b_V = 4*k1_V^3*k3_V - 3*k1_V^2*k2_V^2 - 4*k2_V^3 + k3_V^2 + 6*k1_V*k2_V*k3_V
alpha1_V = (k1_V*k2_V + k3_V - sqrt(b_V))/(2*(k1_V^2 + k2_V))
alpha2_V = (k1_V*k2_V + k3_V + sqrt(b_V))/(2*(k1_V^2 + k2_V))
%beta_V = (2*k1_V^3 + 3*k1_V*k2_V + k3_V - sqrt(b_V))/(sqrt(b_V))
beta_V = (k1_V+alpha2_V)/(alpha1_V-alpha2_V)

T1_V = -(t1_Vc - 0.01)/log(alpha1_V)
T2_V = -(t1_Vc - 0.01)/log(alpha2_V)
T3_V = (beta_V*(T1_V - T2_V)) + T1_V

G_V = tf(K_Vc*[T3_V 1], conv([T1_V 1],[T2_V 1]) )
%--------------------------------------------------------------------------
%RESOLVIENDO (EC23) PARA I
b_I = 4*k1_I^3*k3_I - 3*k1_I^2*k2_I^2 - 4*k2_I^3 + k3_I^2 + 6*k1_I*k2_I*k3_I
alpha1_I = (k1_I*k2_I + k3_I - sqrt(b_I))/(2*(k1_I^2 + k2_I))
alpha2_I = (k1_I*k2_I + k3_I + sqrt(b_I))/(2*(k1_I^2 + k2_I))
%beta_I = (2*k1_I^3 + 3*k1_I*k2_I + k3_I - sqrt(b_I))/(sqrt(b_I))
beta_I = (k1_I+alpha2_I)/(alpha1_I-alpha2_I)

T1_I = -(t1_I - 0.01)/log(alpha1_I)
T2_I = -(t1_I - 0.01)/log(alpha2_I)
T3_I = (beta_I*(T1_I - T2_I)) + T1_I

G_I = tf( K_I*[T3_I 1], conv([T1_I 1],[T2_I 1]) )
%--------------------------------------------------------------------------
%OBTENCIÓN DE PARÁMETROS---------------------------------------------------
C = G_I.num{1}(2)
L = (G_I.den{1}(1))/C
R = (G_I.den{1}(2))/C
%--------------------------------------------------------------------------
%COMPARACIÓN DE GRAFOS-----------------------------------------------------
[Vobt, tobt] = lsim(G_V, u, t);
[Iobt, tobt] = lsim(G_I, u, t);
figure(3)
subplot(3,1,1)
plot(t,u);
title('Entrada'); xlabel('Tiempo [s]'); ylabel('Tensión [V]');
subplot(3,1,2)
plot(t,Vc,t,Vobt);
title('Capacitor');xlabel('Tiempo [s]'); ylabel('Tensión [V]');
subplot(3,1,3)
plot(t,I,t,Iobt);
title('Corriente'); xlabel('Tiempo [s]'); ylabel('Corriente [I]');
%--------------------------------------------------------------------------
%MATRICES------------------------------------------------------------------
A = [-R/L -1/L; 1/C 0];
B = [1/L; 0];
C = [1 0];
D = 0;

G1 = ss(A,B,C,D);
[Icomp,t] = lsim(G1,u,t);

figure(4)
plot(t, Icomp, t, I);grid on; hold on; title('Comparación de corriente');
legend({'i(t)RLC calculados','i(t) de excel'},'Location','southeast')

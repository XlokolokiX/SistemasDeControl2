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

%LECTURA EXCEL-------------------------------------------------------------
Datos = xlsread('Curvas_Medidas_Motor_2024.xls',1);
t = Datos(1:end,1);
Wm = Datos(1:end,2);
I = Datos(1:end,3);
Va = Datos(1:end,4);
Tl = Datos(1:end,5);
Ei = 12;
deltaT = Datos(1,1);
delayT = Datos(702,1)
%--------------------------------------------------------------------------
%GR�FICA DEL EXCEL---------------------------------------------------------
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
%--------------------------------------------------------------------------
%M�TODO DE CHEN
%ALGORITMO DE B�SQUEDA DE LOS PUNTOS---------------------------------------
tW_Inicio = 0.03515;
tW_Delta = 0.00005;
tI_Inicio = 0.03515;
tI_Delta = 0.00005;
[~ , Ind1_I] = min(abs(t-tI_Inicio))
[~ , Ind2_I] = min(abs(t-(tI_Inicio + tI_Delta)))
[~ , Ind3_I] = min(abs(t-(tI_Inicio + 2*tI_Delta)))
[~ , Ind1_W] = min(abs(t-tW_Inicio))
[~ , Ind2_W] = min(abs(t-(tW_Inicio + tW_Delta)))
[~ , Ind3_W] = min(abs(t-(tW_Inicio + 2*tW_Delta)))
t1_W = Datos(Ind1_W,1); W1 = Datos(Ind1_W,2);
t2_W = Datos(Ind2_W,1); W2 = Datos(Ind2_W,2);
t3_W = Datos(Ind3_W,1); W3 = Datos(Ind3_W,2);
t1_I = Datos(Ind1_I,1); I1 = Datos(Ind1_I,3);
t2_I = Datos(Ind2_I,1); I2 = Datos(Ind2_I,3);
t3_I = Datos(Ind3_I,1); I3 = Datos(Ind3_I,3);

%Visualizaci�n de los puntos
zoom = [695 715];
figure(2);
subplot(2,1,1)
plot(t(zoom(1):zoom(2)) ,Wm(zoom(1):zoom(2)), [t1_W t2_W t3_W],[W1 ,W2 ,W3],'+'); grid;
ylabel('[rad/s]'); xlabel('[s]'); title('Velocidad Angular'); ylim([0 250]);
subplot(2,1,2)
plot( t(zoom(1):zoom(2)) ,I(zoom(1):zoom(2)) , [t1_I t2_I t3_I],[I1 ,I2 ,I3],'+'); grid;
ylabel('[A]'); xlabel('[s]'); title('Corriente');
%--------------------------------------------------------------------------
%GANANCIAS-----------------------------------------------------------------
k_W = Datos(16683,2)/Ei; k_I = Datos(16572,3)/Ei;
k1_W = (W1/(Ei*k_W))-1; k2_W = (W2/(Ei*k_W))-1; k3_W = (W3/(Ei*k_W))-1;
k1_I = (I1/(Ei*k_I))-1; k2_I = (I2/(Ei*k_I))-1; k3_I = (I3/(Ei*k_I))-1;
%--------------------------------------------------------------------------
%RESOLVIENDO (EC23) PARA OMEGA MOT-----------------------------------------
b_W = 4*(k1_W^3)*k3_W - 3*(k1_W^2)*(k2_W^2) - 4*(k2_W^3) + (k3_W^2) + 6*k1_W*k2_W*k3_W;
alpha1_W = (k1_W*k2_W + k3_W - sqrt(b_W))/(2*(k1_W^2 + k2_W));
alpha2_W = (k1_W*k2_W + k3_W + sqrt(b_W))/(2*(k1_W^2 + k2_W));
%beta_W = (2*k1_W^3 + 3*k1_W*k2_W + k3_W - sqrt(b_W))/(sqrt(b_W));
beta_W = (k1_W + alpha2_W)/(alpha1_W - alpha2_W);

T1_W = -(t1_W - delayT)/log(alpha1_W);
T2_W = -(t1_W - delayT)/log(alpha2_W);
T3_W = beta_W*(T1_W-T2_W) + T1_W;

G_W = tf(k_W*[T3_W 1],conv([T1_W 1],[T2_W 1]))
%--------------------------------------------------------------------------
%RESOLVIENDO (EC23) PARA CORRIENTE-----------------------------------------
b_I = 4*k1_I^3*k3_I - 3*k1_I^2*k2_I^2 - 4*k2_I^2 + k3_I^3 + 6*k1_I*k2_I*k3_I;
alpha1_I = (k1_I*k2_I + k3_I - sqrt(b_I))/(2*(k1_I^2 + k2_I));
alpha2_I = (k1_I*k2_I + k3_I + sqrt(b_I))/(2*(k1_I^2 + k2_I));
%beta_I = (2*k1_I^3 + 3*k1_I*k2_I + k3_I - sqrt(b_I))/(sqrt(b_I));
beta_I = (k1_I + alpha2_I)/(alpha1_I - alpha2_I);

T1_I = -(t1_I - delayT)/log(alpha1_I);
T2_I = -(t1_I - delayT)/log(alpha2_I);
T3_I = beta_I*(T1_I - T2_I) + T1_I;

G_I = tf(k_I*[T3_I 1],conv([T1_I 1],[T2_I 1]))
%--------------------------------------------------------------------------
%COMPARACI�N CON EXCEL (TRANSITORIO)---------------------------------------
[Wobt, tobt] = lsim(G_W, Va(1:1500), t(1:1500));
[Iobt, tobt] = lsim(G_I, Va(1:1500), t(1:1500));

figure(3)
subplot(2,1,1)
plot(tobt,Wobt ,tobt,Wm(1:1500));
title('Velocidad Angular'); xlabel('Tiempo [s]'); ylabel('[Rad/s]'); ylim([0 250]);legend('Obtenida','Excel')
subplot(2,1,2)
plot(tobt,Iobt, tobt,I(1:1500));
title('Corriente'); xlabel('Tiempo [s]'); ylabel('[I]');legend('Obtenida','Excel')
%--------------------------------------------------------------------------
%OBTENCI�N DE PAR�METROS---------------------------------------------------
ki = (G_W.num{1}(3))
J = G_I.num{1}(2)
B = G_I.num{1}(3)
L = G_I.den{1}(1)/J
R = (G_I.den{1}(2)-L*B)/J
km = (1-R*B)/ki
%--------------------------------------------------------------------------


%PAR�METROS DEL MOTOR (RECALCULADO PARA EL TORQUE)-------------------------
L = 5e-4; J = 2.35e-9; ki = 0.0094054;                                         
%--------------------------------------------------------------------------
%COMPROBACION DEL MODELO CON LA ENTRADA TORQUE-----------------------------
tFinal = t(end);

%Variables de estados, x1 = ia; x2 = wr; u = Va, Tl
%Sistema M I S O
A = [-R/L -km/L ; ki/J -B/J]
B = [1/L 0 ; 0 -1/J]
C = [0 1]                                                                 %Salida omega_R
D = [0 0]

Ts = 0.5e-5;

%Creo mi propio vector tiempo , TL y Va
N = floor(tFinal/Ts)
tProp = 0:Ts:N*Ts;
tl = linspace(0,0,N);
for i=1: N+1 
    [~ , ind] = min( abs(t-i*Ts) );
    if (ind+1) < length(Tl)
        tl(i) = ( (Tl(ind) +  Tl(ind+1))/2);
    else
        tl(i) = Tl(ind);
    end  
end

va = linspace(0,0,N);
for i=1: N+1 
    [~ , ind] = min( abs(t-i*Ts) );
    if (ind+1) < length(Va)
        va(i) = ( (Va(ind) +  Va(ind+1))/2);
    else
        va(i) = Va(ind);
    end  
end

%Simulaci�n del modelo con las entradas obtenidas
x = [0 0]';
Wm_(1) = 0;
for i=1 : N
    xp = A*x + B*[va(i) tl(i)]';
    x = x + xp*Ts;
    Wm_(i+1) = C*x;
end
figure(4)
plot(t,Wm ,tProp,Wm_);
title('Velocidad Angular (Con perturbaci�n)'); xlabel('Tiempo [s]'); ylabel('[Rad/s]'); ylim([0 250]);legend('Obtenida','Excel')
%--------------------------------------------------------------------------

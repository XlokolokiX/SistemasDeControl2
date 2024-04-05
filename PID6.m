% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Mouton Laudin, Alfonso
% Tp N° 1 - Caso de estudio 1 - 
%
%   Inciso 6  
%    Implementar un PID en tiempo discreto para que el ángulo del motor permanezca en una 
%    referencia de 1radian sometido al torque descripto en la Fig. 1-3. (Tip: partir de KP=0,1; 
%    Ki=0,01; KD=5).    
%%
clc; clear all; close all;

%PARÁMETROS DEL MOTOR------------------------------------------------------
Laa = 4.99e-4; J = 4.1302e-6; Ra = 19.9758; Bm = 7.4654e-16; Ki = 16.5207; Km = 0.0605;     
Ei = 12;                                                                   
%--------------------------------------------------------------------------
%MODELADO EN ESP. ESTADO---------------------------------------------------
%Variables de estados, x1 = ia; x2 = wr; x3 = theta; u = Va, Tl
%Sistema M I S O
A = [-Ra/Laa -Km/Laa 0 ; Ki/J -Bm/J 0 ; 0 1 0]
B = [1/Laa 0 ; 0 -1/J ; 0 0]
C = [0 1 0]                                                                %Salida omega_R
D = [0 0]

[num, den] = ss2tf(A,B,C,D,1);
G = tf(num, den)
[num1, den1] = ss2tf(A,B,C,D,2);
G1 = tf(num1, den1)
%RESPUESTAS DEL SISTEMA----------------------------------------------------
poles = pole(G);
tR = log(0.95)/real(poles(2));
Ts = tR/2
%--------------------------------------------------------------------------
%LECTURA DEL EXCEL---------------------------------------------------------
excel = xlsread('Curvas_Medidas_Motor_2024.xls');
t = excel(1:end,1);
TL = excel(1:end,5);
tFinal = t(length(t));
%--------------------------------------------------------------------------
%TORQUE DESCRITO-----------------------------------------------------------
N = floor(tFinal/Ts)                                                       %N puntos a simular
tsim = 0:Ts:N*Ts;
tl = linspace(0,0,N);
for i=1: N+1 
    [~ , ind] = min( abs(t-i*Ts) );
    tl(i) = 1000*TL(ind);
end
clear TL;
clear t;
%--------------------------------------------------------------------------
%PID DISCRETO--------------------------------------------------------------
Kp = 10; Ki = 4; Kd = 0.5;
A1=((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts) 
B1=(-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts) 
C1=Kd/Ts
%--------------------------------------------------------------------------
%CONDICIONES INICIALES-----------------------------------------------------
I(1) = 0;
W(1) = 0;
Theta(1) = 0;
x = [I(1) W(1) Theta(1)]';
xop = [0 0 0]';
Wreferencia = 1;
%--------------------------------------------------------------------------
%SIMULACIÓN----------------------------------------------------------------
e = zeros(1, N);
u = zeros(1, N);
for i=1 : N
    e(i) = Wreferencia - Theta(i);
    if i==1
        u = u + A1*e(i);
    end
    if i==2
        u = u + A1*e(i) + B1*e(i-1);
    end
    if i>2
        u = u + A1*e(i) + B1*e(i-1) + C1*e(i-2);
    end
    
    xp = A*x + B*[u(i) tl(i)]';
    x = x + xp*Ts;
    I(i+1) = x(1);
    W(i+1) = C*x;
    Theta(i+1) = x(3);
end
%--------------------------------------------------------------------------
%GRÁFICAS------------------------------------------------------------------
subplot(4,1,1)
plot(tsim, W); grid;
title('Motor'); xlabel('Tiempo [s]'); ylabel('Velocidad angular [Rad/s]'); ylim([-10 10]);
subplot(4,1,2)
plot(tsim, I); grid;
title('Armadura'); xlabel('Tiempo [s]'); ylabel('Corriente [A]'); ylim([0 0.1]); grid;
subplot(4,1,3)
plot(tsim, Theta);
title('Motor');xlabel('Tiempo [s]'); ylabel('Angulo [Rad]'); ylim([0 2]); grid;
subplot(4,1,4)
plot(tsim, tl);
title('Entradas'); xlabel('Tiempo [s]'); ylabel('[m.N.cm][Volts]'); legend('Perturbación','Tensión Armadura')
%--------------------------------------------------------------------------

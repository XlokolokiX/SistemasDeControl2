%% Sistemas de Control II -FCEFyN-UNC 
% Profesor: Dr.Ing. Pucheta, Julian
% Alumno: Mouton Laudin, Alfonso
% Tp N° 2 - ITEM 2 - 
%
% Calcular sistema controlador que haga evolucionar al péndulo en el equilibrio estable. 
%Objetivo de control: partiendo de una condición inicial nula en el desplazamiento y el ángulo en pi, hacer 
%que el carro se desplace a 10 metros evitando las oscilaciones de la masa m, considerando que es una 
%grúa. Una vez que el desplazamiento=10 modificar a m a un valor 10 veces mayor y volver al origen evitando oscilaciones.
%  -Considerar que sólo puede medirse el desplazamiento y el ángulo. 
%  -Especificar el rango posible para el tiempo de muestreo para implementar el sistema en un 
% microcontrolador.  
%  -Determinar el efecto de la nolinealidad en la acción de control, descripta en la Fig. 4, y verificar cuál 
% es el máximo valor admisible de ésa no linealidad.
%%
clc; clear all; close all;

%Parámetros del péndulo
m1  = 0.1; m2 = m1*10; F  = 0.1; l = 1.6; g  = 9.8; M  = 1.5;
%MODELADO EN EE (Contínuo)
%1er Recorrido (0-10):
A_c1 = [0     1               0          0;       %x1 = delta - desplazamiento
        0    -F/M            -m1*g/M     0;       %x2 = delta_p
        0       0               0        1;       %x3 = phi - ángulo
        0   -F/(l*M)     -g*(m1+M)/(l*M) 0]       %x4 = phi_p

%2do Recorrido (10-0):
A_c2 = [0      1             0         0;           
       0    -F/M         -m2*g/M       0;              
       0       0             0         1;               
       0   -F/(l*M)  -g*(m2+M)/(l*M)   0]
   
B_c = [   0   ;
         1/M  ;
          0   ;
       1/(l*M)]
  
C_c = [1 0 0 0;          %Wr
       0 0 1 0]          %Theta
   
D_c = [0; 
      0]
  
G_c1 = ss(A_c1 ,B_c ,C_c ,D_c);        % m1
G_c2 = ss(A_c2 ,B_c ,C_c ,D_c);        % m2
%%
%Cálculo de las dinámicas del sistema:
roots_c1 = real(eig(A_c1))
roots_c2 = real(eig(A_c2))

% Dinámica rápida
tR   = log(0.95)/roots_c1(2) 
% Dinámica lenta
tL   = log(0.05)/roots_c1(3) 
%%
%SISTEMA DISCRETIZADO
%Tm = tR/100                       %Tiempo de Muestreo
Tm = 1e-2;

G_d1 = c2d(G_c1, Tm, 'zoh');
A_d1     = G_d1.a;
B_d1     = G_d1.b;
C_d1     = G_d1.c;
D_d1     = G_d1.d;

G_d2 = c2d(G_c2, Tm, 'zoh');
A_d2     = G_d2.a;
B_d2     = G_d2.b;
C_d2     = G_d2.c;
D_d2     = G_d2.d;
%%
% CONTROLABILIDAD
Mc1 = [B_d1 A_d1*B_d1 A_d1^2*B_d1 A_d1^3*B_d1 A_d1^4*B_d1]; 
Mc2 = [B_d2 A_d2*B_d2 A_d2^2*B_d2 A_d2^3*B_d2 A_d2^4*B_d2]; 
rank(Mc1)
rank(Mc2)

% ALCANZABILIDAD
Ma1 = [B_d1 A_d1*B_d1 A_d1^2*B_d1 A_d1^3*B_d1];
Ma2 = [B_d2 A_d2*B_d2 A_d2^2*B_d2 A_d2^3*B_d2]; 
rank(Ma1)
rank(Ma2)
%%
%SISTEMA AMPLIADO
Aa1 = [A_d1 , zeros(4,1) ; -C_d1(1,:)*A_d1, eye(1)];
Ba1 = [B_d1; -C_d1(1,:)*B_d1];

Aa2 = [A_d2 , zeros(4,1) ; -C_d2(1,:)*A_d2, eye(1)];
Ba2 = [B_d2; -C_d2(1,:)*B_d2];

%%
%LQR CONTROLADORES
d1 = [10 10 100 100 0.001];               % delta, delta_p, phi, phi_p, ei       
Q1 = diag(d1);
R1 = 1;
[KLQR1, ~, ~] = dlqr(Aa1, Ba1, Q1, R1); 
K1  = KLQR1(1:4);     
Ki1 = -KLQR1(5);

d2 = [10 10 100 100 0.001];
Q2 = diag(d2);
R2 = 1;
[KLQR2, ~, ~] = dlqr(Aa2, Ba2, Q2, R2); 
K2  = KLQR2(1:4);     
Ki2 = -KLQR2(5);

%%
%LQR OBSERVADORES

%Sistema DUAL:
Ao1 = A_d1';
Bo1 = C_d1';       
Co1 = B_d1';

Ao2 = A_d2';
Bo2 = C_d2';       
Co2 = B_d2';

do1 = [1 50 500 .1]; 
Qo1 = diag(do1); 
Ro1 = diag([.1 .1]);
do2 = [1 50 500 .1]; 
Qo2 = diag(do2); 
Ro2 = diag([.1 .1]);
      
% Cálculo del Observador
Ko1 = (dlqr(Ao1 ,Bo1 ,Qo1 ,Ro1))';
Ko2 = (dlqr(Ao2 ,Bo2 ,Qo2 ,Ro2))';

%%
%SIMULACIÓN
Tf = 20;                %Tiempo Final de la simulación
%Ti = 1e-4;              %Tiempo de Integración Euler
Ti = Tm/20
N = ceil(Tf/Ti)         %Número de puntos a simular
deathZone = 0.5;        %Zona muerta del actuador
psita = 0;
d_pp = 0;
phi_pp = 0;

t       = 0:Ti:N*Ti;
u       = 0:Ti:N*Ti;
d_p     = 0:Ti:N*Ti;
d       = 0:Ti:N*Ti;
phi_p   = 0:Ti:N*Ti;
phi     = 0:Ti:N*Ti;

%CONDICIONES INICIALES
d(1) = 0;
d_p(1) = 0;
phi(1) = pi;
phi_p(1) = 0;

x = [d(1) d_p(1) phi(1) phi_p(1)]';
xop = [0 0 pi 0]';

x_obs = [d(1) d_p(1) phi(1) phi_p(1)]';

K = K1;
Ki = Ki1;
Ko = Ko1;
A = A_d1;
B = B_d1;
m = m1;
ref = 10;
f = 0;

%Como el tiempo de integracion mucho más bajo que el tiempo de
%muestreo,actualizo el valor del controlador cada periodo de integración K
SAMPLE_T = floor(Tm/Ti);
um = 0;
y = [0 0]';
y_obs = [0 0]';

for i=1 : N
    
    %Actualizo la acción de control por cada período de muestreo
    if(SAMPLE_T == 0)

        y = C_d1*x;
        y_obs = C_d1*(x_obs + xop);

        psita_p = ref - y(1);
        psita = psita + psita_p;
    
        %Ley de control
        um = -K*(x_obs - xop) + Ki*psita;  %Con Observador
        %um = -K*(x - xop) + Ki*psita;       %Sin Observador

        %No linealidad del Actuador
        if(abs(um) < deathZone)
            um = 0;               
        else
            um = sign(um)*(abs(um) - deathZone);
        end
    end

    u(i) = um;
    
    %Sistema no lineal:
    d_pp        = (1/(M+m))*(  u(i) - m*l*phi_pp*cos(phi(i)) + m*l*phi_p(i)^2*sin(phi(i)) - F*d_p(i) );
    phi_pp      = (1/l)*( g*sin(phi(i)) - d_pp*cos(phi(i)) );
    d_p(i+1)    = d_p(i) + d_pp*Ti;
    d(i+1)      = d(i) + d_p(i)*Ti;
    phi_p(i+1)  = phi_p(i) + phi_pp*Ti;
    phi(i+1)    = phi(i) + phi_p(i)*Ti;
    
    if(d(i+1) >= 9.99)
        if(f == 0)
            K = K2;
            Ki = Ki2;
            Ko = Ko2;
            A = A_d2;
            B = B_d2;
            m = m2;
            ref = 0;
            f = 1;
        end
    end
    
    x = [d(i+1) d_p(i+1) phi(i+1) phi_p(i+1)]';

    if(SAMPLE_T==0)
        x_obs = A*x_obs + B*u(i) + Ko*(y - y_obs) + xop;
        SAMPLE_T = floor(Tm/Ti);
    end
    SAMPLE_T = SAMPLE_T-1;
end

%%
%GRÁFICOS

figure(1);
subplot(3,2,1); grid on; hold on;
plot(t,phi_p,'LineWidth',1.5);grid on; title('Velocidad angular \phi_p');

subplot(3,2,2); grid on; hold on;
plot(t,phi,'LineWidth',1.5); title('Ángulo \phi');xlabel('Tiempo');

subplot(3,2,3); grid on; hold on;
plot(t,d,'LineWidth',1.5);title('Posición grúa \delta');xlabel('Tiempo');

subplot(3,2,4); grid on; hold on;
plot(t,d_p,'LineWidth',1.5);title('Velocidad de grúa \delta_p');

subplot(3,1,3); grid on; hold on;
plot(t,u,'LineWidth',1.5);title('Acción de control u');xlabel('Tiempo en Seg.');
 
figure(2);
subplot(2,1,1);grid on; hold on;
plot(phi,phi_p,'LineWidth',1.5);
title('Ángulo vs Velocidad angular');
xlabel('Ángulo');ylabel('Velocidad angular');
 
subplot(2,1,2);grid on; hold on;
plot(d,d_p,'LineWidth',1.5);
title('Distancia vs velocidad');
xlabel('Distancia');ylabel('Velocidad');
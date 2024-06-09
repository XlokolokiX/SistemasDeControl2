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

A_c = [-a   a   0  0;       %x1 = alpha (ángulo con la horizontal)
       0    0   1  0;       %x2 = phi    (ángulo de cabeceo)
       w^2 -w^2 0  0;       %x3 = phi_p
       c    0   0  0]       %x4 = altura
   
B_c =  [0;
        0;
        b*w^2;
        0]

C_c = [0   0   0   1;      %Altura
       0   1   0   0]      %Angulo de cabeceo

D_c = [0]

G = ss(A_c, B_c ,C_c , D_c)
%Determinación de los tiempos de dinámica
roots_c = eig(A_c)
tR = abs(log(0.95)/real(roots_c(2)))

%%
%SISTEMA DISCRETIZADO
Tm = tR/12

G_d = c2d(G, Tm , 'zoh');
A_d = G_d.a
B_d = G_d.b
C_d = G_d.c
D_d = G_d.d
%%
%SISTEMA AMPLIADO
Aa = [A_d , zeros(4,1) ; -C_d(1,:)*A_d, eye(1)];
Ba = [B_d; -C_d(1,:)*B_d];

%%
%Prueba de Controlabilidad - Alcanzabilidad
Mc = [Ba  Aa*Ba  Aa^2*Ba Aa^3*Ba Aa^4*Ba];
rank(Mc)
%%
%CONTROLADOR
coef_A = poly(Aa)

W = [coef_A(5) coef_A(4) coef_A(3) coef_A(2) 1;
     coef_A(4) coef_A(3) coef_A(2) 1         0;
     coef_A(3) coef_A(2) 1         0         0;
     coef_A(2) 1         0         0         0;
     1         0         0         0         0];
 
T = Mc*W;
p1 = -15+15*i; p2 = -15-15*i ;p3 = -0.5+0.5*i ;p4 = -0.5-0.5*i;     %Polos en plano S
p5 = -0.2;%Polo del integrador
%Transformación Bilineal (plano z)
polos_d = exp([p1 p2 p3 p4 p5]*Tm)
alfa_i = poly(polos_d)

K = fliplr([alfa_i(2:end) - coef_A(2:end)])*inv(T)

Ki = -K(end)
K  = K(1:4)

%%
%OBSERVADOR
%SISTEMA DUAL
Ao = A_d';
Bo = C_d';
Co = B_d';

do = [1 1 1 1];
Qo = diag(do);
Ro = diag([1 1]);
Ko = dlqr(Ao, Bo, Qo, Ro)'

%%
%SIMULACIÓN
Tf = 70;            %Tiempo Final para la Simulación
Ti = 0.01;        	%Tiempo de Integración Euler
N  = ceil(Tf/Ti)    %Numero de puntos a simular
deathZone = 0.5;


t       = 0:Ti:N*Ti;
u       = zeros(1, N+1);
alpha   = zeros(1, N+1);
phi     = zeros(1, N+1);
phi_p   = zeros(1, N+1);
h       = zeros(1, N+1);

%CONDICIONES INICIALES
alpha(1) = 0;
phi(1) = 0;
phi_p(1) = 0;
h(1) = 500;
x = [alpha(1) phi(1) phi_p(1) h(1)]';
x_obs = [0 0 0 0]';
ref = 100;
psita = 0;

um = 0;
SAMPLE_T = floor(Tm/Ti);
for i=1: N
    
    if(SAMPLE_T == 0)
        
        y_obs = C_d*x_obs;
        y = C_d*x;
    
        psita_p = ref - y(1);
        psita = psita + psita_p;
    
        %um = -K*x + Ki*psita;         %Sin Observador
        um = -K*x_obs + Ki*psita;     %Con Observador
    
        %No linealidad del Actuador
        if(abs(um)<deathZone)
            um = 0;
        else
            um = sign(um)*(abs(um) - deathZone);
        end
    end
    u(i) = um;
    
    %SISTEMA
    alpha_p = a*(phi(i) - alpha(i));
    phi_pp = -w^2*(phi(i) - alpha(i) -b*u(i));
    h_p = c*alpha(i);
    alpha(i+1) = alpha(i) + Ti*alpha_p;
    phi_p(i+1) = phi_p(i) + Ti*phi_pp;
    phi(i+1) = phi(i) + Ti*phi_p(i);
    h(i+1) = h(i) + Ti*h_p;
    
    x = [alpha(i+1) phi(i+1) phi_p(i+1) h(i+1)]';
    
    if(SAMPLE_T == 0)
       x_obs = A_d*x_obs + B_d*u(i) +  Ko*(y - y_obs);
       SAMPLE_T = floor(Tm/Ti);
    end
    SAMPLE_T = SAMPLE_T-1;
end

%%
%GRÁFICAS
figure;
subplot(3,1,1);
plot(t, h, 'LineWidth',1.5); title('Altura del avión'); xlabel('Tiempo [s]'); ylabel('Altura [m]');
grid on;
subplot(3,1,2);
plot(t, phi, 'LineWidth',1.5); title('Cabeceo'); xlabel('Tiempo [s]'); ylabel('Ángulo de Cabeceo [Rad]');
grid on;
subplot(3,1,3);
plot(t, u, 'LineWidth',1.5); title('Acción de control'); xlabel('Tiempo [s]'); ylabel('u');
grid on;
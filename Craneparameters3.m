%Open Loop Indoor crane anti-sway control


%Assumptions:
%Zero disturbances
%Negligible string elongation
%Zero friction
%Small deviation of theta
%Neglect effect of initial jerk 


clc
clear all
close all

global M m l g D t1 t2 Vmax w T thetamax a_max a n t_total V_average V2 T1 T2 T3 T4 a2 V1; 

% Default l=40, M=26 tons, m=67.3 tons, Vmax=3m/s


% Crane constants
M = 26*10^3;  % Mass of trolley (kg)
m = 67.3*10^3;  % Mass of load (kg)
l = 40;  % Length of hoist rope (m)
D = 120;  % Desired travel distance (m)

y=  0.00261*D^2-0.5739*D+39.78;
D=D*(1+y/100)-0.6;

Vmax= 3; %Maximum trolley speed (m/s)
g = 9.81;   % Gravitational acceleration(m/s^2)
w = sqrt(g/l); %natural frequency (rad/s)
T= 2*pi/w; %period (s)
thetamax= 5*(pi/180); %max sway angle(rad)

%Velocity profile parameters(s)
a_max=thetamax*g/2; %Estimated max a for thetamax=5deg (m/s^2)
t1 = T/4;% pulse width/non-pulse width (s)

%Velocity profile parameters for t2>0(s)

n = floor(Vmax/(2*a_max*t1)); %number of +ve/-ve pulse pairs, take ceiling to reach Vmax
V1 = 2*n*t1*a_max; 
a=a_max;
a2=(Vmax-V1)/(2*t1);


t2 = (D-(4*n+1)*t1*V1-10*(t1^2)*a2-2*V1*T)/Vmax;  % width when trolley travelling at Vmax(s)

if t2<0
    
    n=floor(Vmax/(2*a_max*t1))-1;
    V2=2*(n+1)*t1*a_max;
    a=a_max;
    a2=a_max;
    t2 = (D-(4*n+5)*t1*V2)/V2;

    if t2<0
        a2=D/(2*(4*n+5)*(n+1)*(t1^2));  %Corresponding a for desired distance D 
        t2=0;
    end 
end 

t_total = 2*(n+1)*T + t2; %Total time for travel
V_average = D/t_total; %Average velocity

T1= n*T; 
T2= (n+1)*T; 
T3= (n+1)*T+t2;
T4= (n+2)*T+t2;



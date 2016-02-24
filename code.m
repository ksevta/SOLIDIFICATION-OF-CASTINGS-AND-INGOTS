clear all
clc
%Aluminum
Ps = 2.7        %Density of Aluminum gm/cm3
H = 95          %Heat of Fusion cal/gm
Tm = 660        %Melting Temp (Degree celsius)
Ks = 0.5        %Thermal Conductivity (cal/sec)/(cm2 C/cm)
Cs = 0.215      %Specific Heat cal/gm C
alpha_s = Ks/(Ps*Cs);

%Sand Mold
Km = 0.00145     %Thermal Conductivity for mold(cal/sec)/(cm2 C/cm)
Pm = 1.5         %Density of Aluminum gm/cm3
Cm = 0.27        %Specific Heat cal/gm C
alpha_m = Km/(Pm*Cm);


Tm = 1000; % melting Temperature
T0 = 300;  % Ambient Temperature
tn = 1000;  % Total Time  
dt = 1;   %Time Step
xn = 100;  %length of mold
alpha = 1; 
dx = 2; 
k_m = alpha_m*dt/(dx^2);% must be less then or equal to 0.5
k_m = 0.3;
k_s = alpha_s*dt/(dx^2);    % must be less then or equal to 0.5
k_s = 0.2;
T = zeros(xn+1,tn);
T1 = zeros(xn+1,tn);
k = 0.3; 
a = 1:xn-1;

%temperature drop in solid
T(1,1) = T0;
for i = 2:tn
    T(1,i) = T(1,i-1)+2*100/tn;
end
for i = 2 : xn+1
    T(i,1) = Tm;
end

for t = 1:tn-1
 for x = 2:xn-1
         T(x,t+1) = T(x,t) + k_s*(T(x+1,t)-2*T(x,t)+T(x-1,t));
 end 
  T(xn,t+1) = 1000;         %T(xn,t)+k_s*(T0 - 2*T(xn,t)+T(xn-1,t)); 
end
dtx = 2:tn-1
for t = 2:tn-1
 for x = 2:xn-1
    if(T(x,t) > 999)
        dtx(t) = (T(x-1,t)-T(x-2,t))/dt;
        break
    end
 end  
end

xa = (-xn+1:1:xn)

%mold ----------------------       
for i = 1:tn
    T1(1,i) = T(2,i);
end
for i = 2 : xn+1
    T1(i,1) = T0;
end

for t = 1:tn-1
 for x = 2:xn-1
         T1(x,t+1) = T1(x,t) + k_m*(T1(x+1,t)-2*T1(x,t)+T1(x-1,t));
 end
  T1(xn,t+1) = T1(xn,t)+k_m*(T0 - 2*T1(xn,t)+T1(xn-1,t)); 
end

 for j = 1 : dt : tn
    for i = 1:xn   
            a(i) = T1(xn-i+1,j);
            if(i < xn)
                a(i+xn) = T(i+2,j);
            end
    end
    a(2*xn-1) = 1000;
    a(2*xn) = 1000;
    figure(1), clf
    plot(xa,a,'-r+','Linewidth',1,'Markersize',5);
    xlabel('Distance(cm)');
    ylabel('Temperature(degree Celsius)');
    title(['Temperature Vs Distance Curve in Mold at time : ',num2str(j)]);
    drawnow
 end
 l = 1:tn-1
 dst = 2:tn-1     
 for i = 2:tn-1
     dst(i) = dtx(i)*Ks/(H*Ps);
 end
 s = 1:tn-1;
 s(1) = dst(1);
 for i = 2:tn-1
     s(i) = s(i-1)+dst(i);
 end
figure(2), clf
plot(l,s);
xlabel('Time(sec)');
ylabel('Thickness(cm)');
title('Thickness Vs Time');    

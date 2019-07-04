clear all; close all; clc;

% Expert system, Fuzzy Control
% For arbitrary & rational & proper & stable 2nd order system.
T = 0.001; NT=3;
w=logspace(-3,2,1000);
t = 0:T:NT;
N = length(t);
y = zeros(1,N); u = zeros(1,N);
e=zeros(1,N); yd = zeros(1,N);
edot=zeros(1,N);

y_f = zeros(1,N); u_f = zeros(1,N);
e_f=zeros(1,N); yd_f = zeros(1,N);
edot_f=zeros(1,N);

%% Trajectory design
% Design of a desired output
ydmag = [15 10 7 5 3] * pi/180;
ydw = 2*pi*[0.012 0.61 1.32 3.51 6.71];
for k=3:N;
    yd(k) = ydmag(1)*cos(ydw(1)*(k-2)*T)+ydmag(2)*cos(ydw(2)*(k-2)*T)+ydmag(3)*cos(ydw(3)*(k-2)*T)+...
            ydmag(4)*cos(ydw(4)*(k-2)*T)+ydmag(5)*cos(ydw(5)*(k-2)*T) - sum(ydmag);    % A point-to-point trajectory
end
% figure('position',[110 420 400 300],'color','w');
% plot(t,yd*180/pi,'r'); ylabel('\theta (deg)'); xlabel('time (s)')

%% Transfer function (Continuous domain)
s = tf('s');
zeta = 0.5; wn = 15;
G = wn^2/(s^2+2*zeta*wn*s+wn^2)
Gz = c2d(G,T,'zoh');
[a, b] = tfdata(Gz,'v');
% roots_G=pole(G)
% figure('color','w'); step(G)
% figure('color','w'); rlocus(G)
% figure('color','w'); bode(G,w)
% BW_G=bandwidth(G)

%% PD Controller  (Discrete domain)
Kp = 100; Kd = 5;
C=tf([Kd Kp],[1]);
Cz = c2d(C,T,'tustin');
[aC, bC] = tfdata(Cz,'v');
Gc = feedback(C*G,1);
Gcz = feedback(Cz*Gz,1);
Gcz = minreal(Gcz);

%% PD Control Feedback Loop

for k=3:N
    y(k) = -b(2)*y(k-1)-b(3)*y(k-2)+a(1)*u(k)+a(2)*u(k-1)+a(3)*u(k-2);
    e(k) = yd(k) - y(k);
    edot(k) = (e(k)-e(k-1))/T;
    u(k) = Kp*e(k) + Kd*edot(k);
    %     u(k) = -bC(2)*u(k-1)+aC(1)*e(k)+aC(2)*e(k-1);
    
    % Saturation of control input
    if u(k) > 10, u(k) = 10; elseif u(k) < -10, u(k) = -10; end
end

%% [PD Control] Figure
figure('position',[610 350 600 500],'color','w');

subplot(311)
plot(t,yd*180/pi,'r','linewidth',3); hold on; plot(t,y*180/pi,'b','linewidth',2); title('PD control result');
xlabel('Time (sec.)'); ylabel('\theta (deg)');
legend('\theta_d','\theta_{act}')
grid on

subplot(312)
plot(t,e*180/pi,'b','linewidth',2); hold on;
xlabel('Time (sec.)'); ylabel('\theta (deg)');
legend('e')
grid on

subplot(313)
plot(t,u,'b','linewidth',2);
xlabel('Time (sec.)'); ylabel('Voltage(V)');
legend('Control input'); grid on
axis([0 t(end) -5 5])

%% Fuzzy Control Membership function

% $$ Error
e_mbs  = linspace(-3,3,N) * pi/180;
e_mbs0 = linspace(-0.5,0.5,5) * pi/180;  % NH, NL, ZE, PL, PH

% $$ Error rate (Error dot)
edot_mbs  = linspace(-20,20,N) * pi/180;
edot_mbs0 = linspace(-15,15,3) * pi/180;  % N, Z, P
% The lower distance between mbs0, the higher sensitivity (Same as D gain ก่)

% $$ Output (u)
u_mbs  = linspace(-3,3,N);
u_mbs0 = linspace(-2.5,2.5,6);  % NH, NL, ZN, ZP, PL, PH

% Membership function
global Point; Point = [2.1 0] * pi/180;
[e_mbsfn, e_sum, se] = Gauss_mbs(e_mbs,e_mbs0,N);
[edot_mbsfn, edot_sum, sedot] = Gauss_mbs(edot_mbs,edot_mbs0,N);
[u_mbsfn, u_sum, su] = Gauss_mbs(u_mbs,u_mbs0,N);

[Point_e_mbsfn, Point_e_sum, se] = Gauss_mbs(Point(1),e_mbs0,1);
[Point_edot_mbsfn, Point_edot_sum, sedot] = Gauss_mbs(Point(2),edot_mbs0,1);

% Fuzzy Rule
[num,sum_num,mu] = fuzzyrule_specific(Point_e_mbsfn,Point_edot_mbsfn,u_mbs0);
% [num,sum_num,mu] = fuzzyrule(Point_e_mbsfn,Point_edot_mbsfn,u_mbs0);
[Point_u_mbsfn, Point_u_sum, su] = Gauss_mbs(mu,u_mbs0,1);

%% [Membership function] Figure
figure('color','w')

subplot(311);
for k=1:length(e_mbs0), plot(e_mbs*180/pi,e_mbsfn(:,k),'linewidth',2); hold on; end
for k=1:length(e_mbs0), plot(Point(1)*180/pi,Point_e_mbsfn(:,k),'ro','linewidth',2); hold on; end
title(sprintf('e_{mbs}: [%1.4f  %1.4f  %1.4f  %1.4f  %1.4f]',Point_e_mbsfn));
ylabel('\mu_{In1}(x)'); xlabel('Error (deg)')

subplot(312);
for k=1:length(edot_mbs0), plot(edot_mbs*180/pi,edot_mbsfn(:,k),'linewidth',2); hold on; end
for k=1:length(edot_mbs0), plot(Point(2)*180/pi,Point_edot_mbsfn(:,k),'ro','linewidth',2); hold on; end
title(sprintf('edot_{mbs}: [%1.4f  %1.4f  %1.4f  %1.4f  %1.4f]',Point_edot_mbsfn));
ylabel('\mu_{In2}(x)'); xlabel('Error rate (deg/s)')

subplot(313);
for k=1:length(u_mbs0), plot(u_mbs,u_mbsfn(:,k),'linewidth',2); hold on; end
for k=1:length(u_mbs0), plot(mu,Point_u_mbsfn(:,k),'ro','linewidth',2); hold on; end
title(sprintf('u = %2.2f',mu));
ylabel('\mu_{Out}(x)'); xlabel('u (V)')

drawnow;
%% Fuzzy Control feedback loop
for k=3:N
    y_f(k) = -b(2)*y_f(k-1)-b(3)*y_f(k-2)+a(1)*u_f(k)+a(2)*u_f(k-1)+a(3)*u_f(k-2);
    e_f(k) = yd(k) - y_f(k);
    edot_f(k) = (e_f(k)-e_f(k-1))/T;
    
    % Membership function
    [Point_e_mbsfn, Point_e_sum] = Gauss_mbs_sinput(e_f(k),e_mbs0,1,se);
    [Point_edot_mbsfn, Point_edot_sum] = Gauss_mbs_sinput(edot_f(k),edot_mbs0,1,sedot);
    
    % Funzzy Rule
        [num,sum_num,mu] = fuzzyrule_specific(Point_e_mbsfn,Point_edot_mbsfn,u_mbs0);
%     [num,sum_num,mu] = fuzzyrule(Point_e_mbsfn,Point_edot_mbsfn,u_mbs0);
    if sum_num==0, mu = 0; end
    u_f(k) = mu;
    
    % Saturation of control input
    if u_f(k) > 10, u_f(k) = 10; elseif u_f(k) < -10, u_f(k) = -10; end
%     k,Point_edot_mbsfn
end

%% [Fuzzy Control feedback loop] Figure
figure('position',[1210 350 600 500],'color','w');

subplot(311)
plot(t,yd*180/pi,'r','linewidth',3); hold on; plot(t,y_f*180/pi,'b','linewidth',2); title('Fuzzy result');
xlabel('Time (sec.)'); ylabel('\theta (deg)');
legend('\theta_d','\theta_{act}')
grid on

subplot(312)
plot(t,e_f*180/pi,'b','linewidth',2); hold on; plot(t,edot_f*180/pi,'r','linewidth',2);
xlabel('Time (sec.)'); ylabel('\theta (deg)');
legend('e')
grid on

subplot(313)
plot(t,u_f,'b','linewidth',2);
xlabel('Time (sec.)'); ylabel('Voltage(V)');
legend('Control input'); grid on
axis([0 t(end) -5 5])


%% [Combined feedback loop] Figure
figure('position',[1210 350 600 500],'color','w');

subplot(311)
plot(t,yd*180/pi,'g','linewidth',4); hold on; 
plot(t,y*180/pi,'b','linewidth',2); hold on;
plot(t,y_f*180/pi,'r-.','linewidth',1);
title('Ordinary');
xlabel('Time (sec.)'); ylabel('\theta (deg)');
legend('\theta_d','\theta_{act PD}','\theta_{act Fuzzy}')
grid on

subplot(312)
plot(t,e*180/pi,'b','linewidth',2); hold on; 
plot(t,e_f*180/pi,'r-.','linewidth',2);
xlabel('Time (sec.)'); ylabel('\theta (deg)');
legend('PD','Fuzzy')
grid on

subplot(313)
plot(t,u,'b','linewidth',3); hold on;
plot(t,u_f,'r-.','linewidth',1);
xlabel('Time (sec.)'); ylabel('Voltage(V)');
legend('u_{PD}','u_{Fuzzy}'); grid on
axis([0 t(end) -5 5])
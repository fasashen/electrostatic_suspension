clear all
clc
close('all')

eps0 = 8.8541878e-12;

a = 0.05;
b = 1;
h = 30E-06;
S = a*b;
g = 9.8;
m = 3;

L = 10e-3;
C = eps0*S/h;
U0 = 10;
R = 10;
f = 1e4;            % frequency in Hz

w = sqrt(1/L/C);    % circular frequency
n = 0.5*R/L;

A = U0/L/sqrt((w^2 - (2*pi*f)^2)^2 + 4*n^2*(2*pi*f)^2)/C;
phi = -atan( R/((2*pi*f)*L-1/w/C)  );

Afun = @(f) U0/L./sqrt((w^2 - (2*pi*f).^2).^2 + 4*n^2*(2*pi*f).^2)/C;
Phifun = @(f) -atan(2*n*2*pi*f./(w^2 - (2*pi*f).^2));
xi = 0.5*R*sqrt(C/L);
wd = w*sqrt(1-xi^2);

disp(['A = ' num2str(A)])
disp(['Phi = ' num2str(phi)])
disp(['wn = ' num2str(w)]);      % natural frequency of undamped system 
disp(['wd = ' num2str(wd)]);     % natural frequency of damped system
disp(['A(wn) = ' num2str(Afun(w))])
disp(['A(wd) = ' num2str(Afun(wd))])

Wv = linspace(0,2*w,300);
gap = linspace(100e-6,300e-6,300);
time = linspace(1.2e-2,1.28e-2,1/f*128/(1/f/16));

% figure(1)
% subplot(3,1,1)
% plot(Wv,Afun(Wv))
% subplot(3,1,2)
% plot(Wv,Phifun(Wv))
% subplot(3,1,3)
% LCRcircuit = @(t,Y) [Y(2); 1/L*(-R*Y(2) + U0*sin(2*pi*f*t) -1/C*Y(1))];
% [t, y] = ode23t(LCRcircuit,[0 64*1/f],[0 0]);
% % t = linspace(0,1/W*20,300);
% % plot(t,A*sin(W*t - phi))
% plot(t,1/C*y(:,1))

gap = linspace(10e-6,60e-6,300);

figure(1)
subplot(3,1,1)
q2_fun = @(h) U0^2/L^2./(((1/L./(eps0*S./h)) - (2*pi*f).^2).^2 + 4*n^2*(2*pi*f).^2)/2;
MotionEquation = @(t,Y) [Y(2); -g+q2_fun(Y(1))^2/(2*eps0*S*m)];
[t, y] = ode23t(MotionEquation,[0 64*3*1/f],[h 0]);
plot(t,y(:,1))
subplot(3,1,2)
plot(gap,q2_fun(gap).*gap/eps0/S)
subplot(3,1,3)
LCRcircuit = @(t,Y) [Y(2); 1/L*(-R*Y(2) + U0*sin(2*pi*f*t) -1/C*Y(1))];
[t, y] = ode23t(LCRcircuit,[0 64*1/f],[0 0]);
plot(t,1/C*y(:,1))


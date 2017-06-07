clear all
clc
close('all')


eps0 = 8.8541878e-12;

a = 0.005;
b = 0.013;
h = a-b;
d = 3e-6;
S = 2*3.1415*a*h;

L = 20e-3;
R = 400;
C = eps0*S/d;

U0 = 10;
f = 40000;            % frequency in Hz

w = sqrt(1/L/C);    % circular frequency
n = 0.5*R/L;
g = 9.8;
m = 3;


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
gap = linspace(7e-6,10e-6,100);
time = linspace(1.2e-2,1.28e-2,1/f*128/(1/f/16));

Cfun = @(y) eps0*S./y;
peak_coef = @(y) (2*pi*f)*L-1./(2*pi*f)./Cfun(y);
disp([num2str((2*pi*f)*L) ' - ' num2str(1./(2*pi*f)./Cfun(d))]);
disp(['Peak = ' num2str(peak_coef(d))]);
Ffun = @(y) (U0./((2*pi*f).*(R^2+(peak_coef(y)).^2).^0.5)).^2./(2*eps0*S)./2;

force = Ffun(gap);


figure(1)

subplot(4,1,1)
LCRcircuit = @(t,Y) [Y(2); 1/L*(-R*Y(2) + U0*sin(2*pi*f*t) -1/C*Y(1))];
[t, y] = ode23t(LCRcircuit,[0 64*1/f],[0 0]);

plot(t,1/C*y(:,1))
subplot(4,1,2)
plot(gap,force)

sys_static = @(t,Y) [Y(2); 1/L*(U0*sin(2*pi*f*t) - Y(1).*d/S/eps0 - R*Y(2))];
[t, y] = ode23t(sys_static,[0 64*1/f],[0 0]);
u = y(:,1)./C;
f_func = y(:,1).^2./(2*eps0*S);
subplot(4,1,3)
plot(t,u)
subplot(4,1,4)
plot(t,f_func)


% gap = linspace(10e-6,60e-6,300);
% 
% figure(1)
% subplot(3,1,1)
% q2_fun = @(h) U0^2/L^2./(((1/L./(eps0*S./h)) - (2*pi*f).^2).^2 + 4*n^2*(2*pi*f).^2)/2;
% MotionEquation = @(t,Y) [Y(2); -g+q2_fun(Y(1))^2/(2*eps0*S*m)];
% [t, y] = ode23t(MotionEquation,[0 64*3*1/f],[h 0]);
% plot(t,y(:,1))
% subplot(3,1,2)
% plot(gap,q2_fun(gap).*gap/eps0/S)
% subplot(3,1,3)
% LCRcircuit = @(t,Y) [Y(2); 1/L*(-R*Y(2) + U0*sin(2*pi*f*t) -1/C*Y(1))];
% [t, y] = ode23t(LCRcircuit,[0 64*1/f],[0 0]);
% plot(t,1/C*y(:,1))


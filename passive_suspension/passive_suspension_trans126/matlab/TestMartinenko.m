clear all
clc
eps0 = 8.854e-12;
L = 20e-4;
R = 300;
u0 = 10;
w = 40000;
S = 0.005;
g = 9.8;
m = 0.122395/2/g;
gap = 3e-6;


sys1 = @(t,Y) [Y(2); (1/L*(u0*sin(2*pi*w*t) - Y(1).*Y(3)./S/eps0 - R.*Y(2)));...
    Y(4); (1/(S*2*eps0*m)*Y(1).^2 - g)];

[t, y] = ode23t(sys1,[0 32*1/w],[0 0 gap 0]);

e1 = y(:,1);
y1 = y(:,3);
subplot(2,1,1)
plot(t,y1)
subplot(2,1,2)
plot(t,e1)


% C = eps0*S/gap;
% 
% w0 = 1/sqrt(L*C);
% sys_static = @(t,Y) [Y(2); 1/L*(u0*sin(2*pi*w*t) - Y(1).*gap/S/eps0 - R*Y(2))];
% [t, y] = ode23t(sys_static,[0 16*1/w],[0 0]);
% u = y(:,1)./C;
% f_func = y(:,1).^2./(2*eps0*S);
% subplot(2,1,1)
% plot(t,u)
% subplot(2,1,2)
% plot(t,f_func)

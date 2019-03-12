% 'Robot Dynamics' HW_3 Problem 2
clc;clear all;close all;
syms m1 m2 l1 l2 a1 a2 t1(t) t2(t) I1 I2 g

%%
%Link 1
x1 = l1*cos(t1);
y1 = l1*sin(t1);
x1d = diff(x1,t);
y1d = diff(y1,t);
v1_sq = x1d^2 + y1d^2;
    %Kinetic energy for link1
k1 = (1/2)*m1*v1_sq + (1/2)*I1*diff(t1(t),t)^2
    %Potential energy for link1
p1 = m1*g*y1

%%
%Link 2
x2 = a1*cos(t1) + l2*cos(t1+t2);
y2 = a1*sin(t1) + l2*sin(t1+t2);
x2d = diff(x2,t);
y2d = diff(y2,t);
v2_sq = x2d^2 + y2d^2;
    %Kinetic energy for link2
k2 = (1/2)*m2*v2_sq + (1/2)*I2*(diff(t1(t),t)+diff(t2(t),t))^2
    %Potential energy for link2
p2 = m2*g*y2

%%
%Total kinetic energy
k = k1 + k2;
k = simplify(k)

%%
%Total potential energy
p = p1 + p2;
p = simplify(p)

%%
%Lagrangian, L
L = k-p

%%
%To find torque vector

syms t1d t2d th1 th2

%Substituting terms for allowing differentiation
L = subs(L,[diff(t1(t),t), diff(t2(t),t)],[t1d, t2d]);
L = subs(L,[t1(t), t2(t)],[th1, th2]);

% Partial derivatives for tau1
partial_wrt_t1d = diff(L,t1d);
partial_wrt_t1 = diff(L,th1);

% Partial derivatives for tau2
partial_wrt_t2d = diff(L,t2d);
partial_wrt_t2 = diff(L,th2);

%Reversing substitution to go back to original form
partial_wrt_t1d = subs(partial_wrt_t1d,[t1d, t2d],[diff(t1(t),t), diff(t2(t),t)]);
partial_wrt_t1d = subs(partial_wrt_t1d,[th1 th2],[t1(t) t2(t)]);

partial_wrt_t1 = subs(partial_wrt_t1,[th1, th2],[t1(t), t2(t)]);
partial_wrt_t1 = subs(partial_wrt_t1,[t1d, t2d],[diff(t1(t),t), diff(t2(t),t)]);

partial_wrt_t2d = subs(partial_wrt_t2d,[t1d, t2d],[diff(t1(t),t), diff(t2(t),t)]);
partial_wrt_t2d = subs(partial_wrt_t2d,[th1, th2],[t1(t), t2(t)]);

partial_wrt_t2 = subs(partial_wrt_t2,[th1, th2],[t1(t), t2(t)]);
partial_wrt_t2 = subs(partial_wrt_t2,[t1d, t2d],[diff(t1(t),t), diff(t2(t),t)]);

% Time derivative for tau1
time_derivative_tau1 = diff(partial_wrt_t1d,t);
% Time derivative for tau2
time_derivative_tau2 = diff(partial_wrt_t2d,t);

%%
%Final equation for tau1
tau1 = time_derivative_tau1 - partial_wrt_t1;
tau1 = simplify(tau1)
%Final equation for tau2
tau2 = time_derivative_tau2 - partial_wrt_t2;
tau2 = simplify(tau2)

%The torque vector
Tau = [tau1; tau2]

%%
%Performing substitutions again, to find inertia, gravity, and coriolis matrices.
syms t1dd t2dd
Tau = subs(Tau,[diff(t1(t),t,t), diff(t2(t),t,t)],[t1dd, t2dd]);
Tau = subs(Tau,[diff(t1(t),t), diff(t2(t),t)],[t1d, t2d]);
disp('Your new tau')
disp(Tau)

%%
%Inertia matrix
IM = equationsToMatrix(Tau,[t1dd t2dd])

%%
%Gravity matrix
GM = equationsToMatrix(Tau,g)

%%
%Coriolis matrix
CM = simplify(expand(Tau - IM*[t1dd;t2dd] - GM*g))
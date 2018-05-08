function [rates, sol, t] = bg_sim(v_e, W, c, u0, T, dt, da)
%cortical - e, D1 - d1, D2 - d2, GPi - p1, GPe - p2, STN - stn
% spikes/s
S_p1 = 250.0;
S_p2 = 300.0;  S_d1 = 65.0;  S_d2 = 65.0;  S_stn = 500.0;
% mV^-1
k_p1 = 0.2; k_p2 = 0.2;  k_d1 = 0.3;  k_d2 = 0.3;  k_stn = 0.2;
% mV
V_p1 = 10.0;  V_p2 = 9.0;  V_d1 = 19.0;  V_d2 = 19.0;  V_stn = 10.0;
% decay and rise time ants of membrane s^-1
a = 160.0;  b = 640.0;
% mV s          p1     p2   d1    d2   stn    e
v =  [[    0  -0.03  -0.1   0   0.3   0  ]
        [    0   -0.1    0  -0.3  0.3   0  ]
        [    0     0     0    0    0   1.0/da ]
        [    0     0     0    0    0   0.7*da ]
        [    0  -0.04    0    0    0   0.1 ]
        [    0     0     0    0    0    0  ]]';

delay =  [[    0    1.0   1.0   0   1.0   0  ]
      [    0     0     0   1.0  1.0   0  ]
      [    0     0     0    0    0   2.0 ]
      [    0     0     0    0    0   2.0 ]
      [    0    1.0    0    0    0   1.0 ]
      [    0     0     0    0    0    0  ]]' * 0.001;

tspan = [0.0 T];
[t, sol] = ode1(@bg_model,tspan,u0,dt);
rates = [sigma(sol(:,1), V_p1, k_p1, S_p1) sigma(sol(:,3), V_p2, k_p2, S_p2) sigma(sol(:,5), V_d1,  k_d1,  S_d1) sigma(sol(:,7), V_d2,  k_d2,  S_d2) sigma(sol(:,9), V_stn, k_stn, S_stn)];

function S = sigma(v, V, k, S_max)
  S = S_max./(1+exp(k*(V - v)));
end

function S = dsigma(v, V, k, S_max)
  S = S_max*k*exp(k*(V - v))./(1+exp(k*(V - v))).^2;
end

function du = bg_model(t,u)
    du = zeros(13,1);
    dW = -c*(u(11) - u0(11)) + sqrt(dt)*W(t)/dt;
    %First order approximation for delays
    du(1) = u(2);
    du(2) = (a*b)*v(5,1)*(sigma(u(9), V_stn, k_stn, S_stn) - delay(5,1)*u(10)*dsigma(u(9), V_stn, k_stn, S_stn)) + ...
            (a*b)*v(2,1)*(sigma(u(3), V_p2, k_p2, S_p2) - delay(2,1)*u(4)*dsigma(u(3), V_p2, k_p2, S_p2)) + ...
            (a*b)*v(3,1)*(sigma(u(5), V_d1,  k_d1,  S_d1) - delay(3,1)*u(6)*dsigma(u(5), V_d1,  k_d1,  S_d1)) - ...
            (b + a)*u(2) - (a*b)*u(1);
    du(3) = u(4);
    du(4) = (a*b)*v(5,2)*(sigma(u(9), V_stn, k_stn, S_stn) - delay(5,2)*u(10)*dsigma(u(9), V_stn, k_stn, S_stn)) + ...
            (a*b)*v(4,2)*(sigma(u(7), V_d2,  k_d2,  S_d2) - delay(4,2)*u(8)*dsigma(u(7), V_d2,  k_d2,  S_d2)) + ...
            (a*b)*v(2,2)*(sigma(u(3), V_p2, k_p2, S_p2) - delay(2,2)*u(4)*dsigma(u(3), V_p2, k_p2, S_p2)) - ...
            (b + a)*u(4) - (a*b)*u(3);
    du(5) = u(6);
    du(6) = (a*b)*v(6,3)*u(11) - (b + a)*u(6) - (a*b)*u(5);
    du(7) = u(8);
    du(8) = (a*b)*v(6,4)*u(12) - (b + a)*u(8) - (a*b)*u(7);
    du(9) = u(10);
    du(10) = (a*b)*v(6,5)*u(13) + ...
             (a*b)*v(2,5)*(sigma(u(3), V_p2, k_p2, S_p2) - delay(2,5)*u(4)*dsigma(u(3), V_p2, k_p2, S_p2)) - ...
             (b + a)*u(10) - (a*b)*u(9);
    du(11) = v_e(t-delay(6,3)) + dW;
    du(12) = v_e(t-delay(6,4)) + dW;
    du(13) = v_e(t-delay(6,5)) + dW;
end
end
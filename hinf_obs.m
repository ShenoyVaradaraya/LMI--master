clc;clear;
A = [-1.3410  0.9933 0  -0.1689 -0.2518;...
     43.2230 -0.8693 0 -17.2510 -1.5766;...
      1.3410  0.0067 0   0.1689  0.2518;...
         0      0    0 -20.0000   0;...
         0      0    0    0     -20.0000];
B1 = [0 0; 0 0; 0 0; 20 0; 0 20];  
C2 = [1 0 0 0 0;...
     0 0 1 0 0];
C1 = [1 0; 0 0; 0 0 ;0 0; 0 0];
D1 = [1 0;0 0];
D2 = 0;
P = sdpvar(5);
W = sdpvar(5,2);
gamma = sdpvar(1);

hinf_LMI = [[(P*A + A'*P - C2'*W'-W*C2) (P*B1-W*D2) C1;
              (P*B1-W*D2)' -gamma*eye(2) D1;
              C1' D1 -gamma*eye(2)] <= 0];
hinf_LMI = [hinf_LMI, P >= 0];
opt = sdpsettings('verbose',0,'solver','mosek');
optimize(hinf_LMI,gamma,opt);
P = value(P);
W = value(W);
L = inv(P)*W
fprintf('\n Gamma : %f\n',value(gamma))
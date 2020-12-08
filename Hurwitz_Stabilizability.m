clc;clear;
A =  [ -1 0 1 1 0 0 -1 2;
    0 0 1 0 -3 -1 0 0 ;
    -1 -2 1 0 -3 -1 0 0;
    0 0 0 0 1 0 0 0;
    0 0 -1 0 -3 0 4 2;
    0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 1;
    -1 2 0 0 1 -1 -1 2;];
B = [ 1 zeros(1,3);zeros(1,4); 0 1 zeros(1,2); zeros(1,4); zeros(1,2) 1 0; zeros(3,4)];
C = [ 17 23 3 14 -15 5 25 5; 3 -3 2 3 -2 2 -3 -6; -2 1 0 -1 2 0 1 4];

P = sdpvar(8);
W = sdpvar(4,8);
Hurwitz_LMI = [A*P + P*A' + B*W + W'*B' <= 0];
opt = sdpsettings('verbose',0,'solver','mosek');
optimize(Hurwitz_LMI,[],opt);
P = value(P);
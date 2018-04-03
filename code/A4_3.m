close all
clear 

R1 = 1;
G1 = 1/R1;
c = 0.25;
R2 = 2;
G2 = 1/R2;
L = 0.2;
R3 = 10;
G3 = 1/R3;
a = 100;
R4 = 0.1;
G4 = 1/R4;
Ro = 1000;
Go = 1/Ro;
Vin = 1;

C = [0 0 0 0 0 0 0;
    -c c 0 0 0 0 0;
    0 0 -L 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;];

G = [1 0 0 0 0 0 0;
    -G2 G1+G2 -1 0 0 0 0;
    0 1 0 -1 0 0 0;
    0 0 -1 G3 0 0 0;
    0 0 0 0 -a 1 0;
    0 0 0 G3 -1 0 0;
    0 0 0 0 0 -G4 G4+Go];

F = [Vin;
    0;
    0;
    0;
    0;
    0;
    0;];

Foff = [Vin-Vin;
        0;
        0;
        0;
        0;
        0;
        0;];

ts = 1000;              % Time step
    
V1 = zeros(7, ts);
Vstart = zeros(7, 1);
dt = 1e-3;
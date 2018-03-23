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
    
V1 = zeros(7, 1000);
Vstart = zeros(7, 1);
dt = 1e-3;
for i = 1:1000

    if i < 30
        V1(:,i) = (C./dt+G)\(Foff+C*Vstart/dt);
    elseif i == 30
        V1(:,i) = (C./dt+G)\(F+C*Vstart/dt);
    else
        V1(:,i) = (C./dt+G)\(F+C*Vpast/dt);
    end
    Vpast = V1(:, i);
    
end
figure(1)
plot(1:1000, V1(7,:))
hold on
plot(1:1000, V1(1,:))
title('Vin and Vout with a step at 0.03s')
xlabel('Time (ms)')
ylabel('Voltage (V)')

V2 = zeros(7, 1000);
Fsin = zeros(7,1);
for j = 1:1000

    Vsin = sin(2*pi*(1/0.03)*j/1000);
    Fsin(1,1) = Vsin;
    if j == 1
        V2(:,j) = (C./dt+G)\(Fsin+C*Vstart/dt);
    else
        V2(:,j) = (C./dt+G)\(Fsin+C*Vpast/dt);
    end
    Vpast = V2(:, j);
        
end
figure(2)
plot(1:1000, V2(7,:))
hold on
plot(1:1000, V2(1,:))
title('Vin and Vout with sin function of f = 1/0.03')
xlabel('Time (ms)')
ylabel('Voltage (V)')

V3 = zeros(7, 1000);
Fguass = zeros(7,1);
for k = 1:1000

    Vgauss = exp(-1/2*((k/1000-0.06)/(0.03))^2);
    Fgauss(1,1) = Vgauss;
    if k == 1
        V3(:,k) = (C./dt+G)\(Fgauss+C*Vstart/dt);
    else
        V3(:,k) = (C./dt+G)\(Fgauss+C*Vpast/dt);
    end
    Vpast = V3(:, k);
        
end
figure(3)
plot(1:1000, V3(7,:))
hold on
plot(1:1000, V3(1,:))
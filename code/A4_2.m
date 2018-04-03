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
for i = 1:ts

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
plot(1:ts, V1(7,:), 'r')
hold on
plot(1:ts, V1(1,:), 'b')
title('Vin and Vout with a step at 0.03s')
xlabel('Time (ms)')
ylabel('Voltage (V)')
grid on

V2 = zeros(7, ts);
Fsin = zeros(7,1);
for j = 1:ts

    Vsin = sin(2*pi*(1/0.03)*j/ts);
    Fsin(1,1) = Vsin;
    if j == 1
        V2(:,j) = (C./dt+G)\(Fsin+C*Vstart/dt);
    else
        V2(:,j) = (C./dt+G)\(Fsin+C*Vpast/dt);
    end
    Vpast = V2(:, j);
        
end
figure(2)
plot(1:ts, V2(7,:), 'r')
hold on
plot(1:ts, V2(1,:), 'b')
title('Vin and Vout with sin function of f = 1/0.03')
xlabel('Time (ms)')
ylabel('Voltage (V)')
grid on

V3 = zeros(7, ts);
Fgauss = zeros(7,1);
for k = 1:ts

    Vgauss = exp(-1/2*((k/ts-0.06)/(0.03))^2);
    Fgauss(1,1) = Vgauss;
    if k == 1
        V3(:,k) = (C./dt+G)\(Fgauss+C*Vstart/dt);
    else
        V3(:,k) = (C./dt+G)\(Fgauss+C*Vpast/dt);
    end
    Vpast = V3(:, k);
        
end
figure(3)
plot(0:ts-1, V3(7,:), 'r')
hold on
plot(0:ts-1, V3(1,:), 'b')
title('Vin and Vout with gaussian function')
xlabel('Time (ms)')
ylabel('Voltage (V)')
grid on

f = (-ts/2:ts/2-1);               % Frequency range

fV1in = fft(V1(1, :));
fV1out = fft(V1(7, :));
fsV1in = fftshift(fV1in);
fsV1out = fftshift(fV1out);
figure(4)
plot(f, abs(fsV1in), 'r')
hold on
plot(f, abs(fsV1out), 'b')
title('Vin and Vout in frequency domain with a step function')
xlabel('frequency (1/ms)')
ylabel('Voltage (V)')
grid on

fV2 = fft(V2.');
fsV2 = fftshift(fV2);
figure(5)
plot(f, abs(fsV2(:, 1)), 'r')
hold on
plot(f, abs(fsV2(:, 7)), 'b')
title('Vin and Vout in frequency domain with a sin function')
xlabel('frequency (1/ms)')
ylabel('Voltage (V)')
grid on

fV3 = fft(V3.');
fsV3 = fftshift(fV3);
figure(6)
plot(f, abs(fsV3(:, 1)), 'r')
hold on
plot(f, abs(fsV3(:, 7)), 'b')
title('Vin and Vout in frequency domain with a gaussian function')
xlabel('frequency (1/ms)')
ylabel('Voltage (V)')
grid on


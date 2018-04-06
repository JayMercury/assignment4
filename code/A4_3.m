%% Part 3
% In part 3, we are adding a current source and a capacitor into the
% circuit and simulate it in order to model thermal noise generated in the resistor R3 and
% bandwidth limit to noise respectively.

% Reset Everything
close all
clear 

% Define Constant
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
cn1 = 1e-5;                 % Default capacitances
cn2 = 1e-12;                % Two new capacitances
cn3 = 5.2e-5;

% Creating Matrices
C1 = [0 0 0 0 0 0 0;
    -c c 0 0 0 0 0;
    0 0 -L 0 0 0 0;
    0 0 0 -cn1 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 -cn1 0 0 0;
    0 0 0 0 0 0 0;];

C2 = [0 0 0 0 0 0 0;
    -c c 0 0 0 0 0;
    0 0 -L 0 0 0 0;
    0 0 0 -cn2 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 -cn2 0 0 0;
    0 0 0 0 0 0 0;];

C3 = [0 0 0 0 0 0 0;
    -c c 0 0 0 0 0;
    0 0 -L 0 0 0 0;
    0 0 0 -cn3 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 -cn3 0 0 0;
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

ts1 = 1000;              % Default Time step
ts2 = 1.9898e4;          % New Time step

V1 = zeros(7, ts1);
Vstart = zeros(7, 1);
dt1 = 1e-3;
dt2 = 1.9898e-4;

% Circuit with Noise simulation with default time step
% Time domain simulation
V1 = zeros(7, ts1);
Fgauss = zeros(7,1);
for i = 1:ts1
    
    Fgauss(1,1) = exp(-1/2*((i/ts1-0.06)/(0.03))^2);
    Fgauss(4,1) = 0.001*randn();
    Fgauss(7,1) = 0.001*randn();
    if i == 1
        V1(:,i) = (C1./dt1+G)\(Fgauss+C1*Vstart/dt1);
    else
        V1(:,i) = (C1./dt1+G)\(Fgauss+C1*Vpast/dt1);
    end
    Vpast = V1(:, i);
        
end
figure(1)
plot(1:ts1, V1(7,:), 'r')
hold on
plot(1:ts1, V1(1,:), 'b')
title('Vin with thermal noise and Vout using gaussian excitation')
xlabel('Time (ms)')
ylabel('Voltage (V)')
grid on

% Frequency domain simulation
f = (-ts1/2:ts1/2-1);               % Frequency range

fV1 = fft(V1.');
fsV1 = fftshift(fV1);
figure(2)
plot(f, abs(fsV1(:, 1)), 'r')
hold on
plot(f, abs(fsV1(:, 7)), 'b')
title('Vin with thermal noise and Vout in frequency domain using gaussian excitation')
xlabel('frequency (1/ms)')
ylabel('Voltage (V)')
grid on

% Time domain simulation with smaller capacitance than default one
V2 = zeros(7, ts1);
Fgauss = zeros(7,1);
for j = 1:ts1
    
    Fgauss(1,1) = exp(-1/2*((j/ts1-0.06)/(0.03))^2);
    Fgauss(4,1) = 0.001*randn();
    Fgauss(7,1) = 0.001*randn();
    if j == 1
        V2(:,j) = (C2./dt1+G)\(Fgauss+C2*Vstart/dt1);
    else
        V2(:,j) = (C2./dt1+G)\(Fgauss+C2*Vpast/dt1);
    end
    Vpast = V2(:, j);
        
end
figure(3)
plot(1:ts1, V2(7,:), 'r')
hold on
plot(1:ts1, V2(1,:), 'b')
title('Vin with smaller Cn and Vout using gaussian excitation')
xlabel('Time (ms)')
ylabel('Voltage (V)')
grid on

% Time domain simulation with slighter bigger capacitance than default one
V3 = zeros(7, ts1);
Fgauss = zeros(7,1);
for k = 1:ts1
    
    Fgauss(1,1) = exp(-1/2*((k/ts1-0.06)/(0.03))^2);
    Fgauss(4,1) = 0.001*randn();
    Fgauss(7,1) = 0.001*randn();
    if k == 1
        V3(:,k) = (C3./dt1+G)\(Fgauss+C3*Vstart/dt1);
    else
        V3(:,k) = (C3./dt1+G)\(Fgauss+C3*Vpast/dt1);
    end
    Vpast = V3(:, k);
        
end
figure(4)
plot(1:ts1, V3(7,:), 'r')
hold on
plot(1:ts1, V3(1,:), 'b')
title('Vin with slightly bigger Cn and Vout using gaussian excitation')
xlabel('Time (ms)')
ylabel('Voltage (V)')
grid on

%% Discussion
% By increasing the capacitance just slightly, the simulation breaks down
% since the circuit is trapped in a feedback loop but nothing is changed
% when the capacitance is too small.

% Time domain simulation with slight smaller time step than default
V4 = zeros(7, ts2);
Fgauss = zeros(7,1);
for m = 1:ts2
    
    Fgauss(1,1) = exp(-1/2*((m/ts2-0.06)/(0.03))^2);
    Fgauss(4,1) = 0.001*randn();
    Fgauss(7,1) = 0.001*randn();
    if m == 1
        V4(:,m) = (C1./dt2+G)\(Fgauss+C1*Vstart/dt2);
    else
        V4(:,m) = (C1./dt2+G)\(Fgauss+C1*Vpast/dt2);
    end
    Vpast = V4(:, m);
        
end
figure(5)
plot(1:ts2, V4(7,:), 'r')
hold on
plot(1:ts2, V4(1,:), 'b')
title('Vin with thermal noise and Vout with smaller timestep')
xlabel('Time (ps)')
ylabel('Voltage (V)')
grid on

%% Discussion
% By just decreasing the time step just by a little bit, the simulation
% also seems to break down since the circuit seems to be trapped in a
% feedback loop

%% Part 4
% In part 4, we would change the voltage source on the output stage from
% the equation $V = a*I3$ to $V = alpha*I3 + beta*I3^2 + gamma*I3^3$, this
% simple change of equation doesn't mean we can just change the matrices we
% have to fit the new equation. We would have to create new matrices for
% the Jacobian method on simulating the new equation in the circuit. This
% would also mean that new equations are iterate and increase the size of
% the matrices. The beta and gamma can be defined by us, but both of them
% has to very large in order have effect on the simulation since if I3 is
% smaller than 1, the square and the cube of I3 will only make it smaller.

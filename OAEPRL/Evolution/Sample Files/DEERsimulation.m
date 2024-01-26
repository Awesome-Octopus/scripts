%%  This simulation shows the dependence of the postions of the spins that feel the
% effect of the pump pulse as a function of the electron-electron coupling and the position of the pump pulse

clear all

%% Sets up the pump pulse timing and electron-electron coupling.
% Changing the position moves the pump pulse
% eecoupling is the dipolar coupling for one orientation

position = 166.67;  % Position of the pump pump with respect to the first echo in ns
eecoupling = 3;  %  Electron-Electron coupling in MHz

%%  Carries out the DEER Pulse sequence

timing = position/1e9;
offset = eecoupling*1e6;
tau = 200e-9;
T = 1e-6;

% Sets up the spins along the Z axis
for i=1:250;
	Spins.degrees(i,:) = [0,0,1];
end

Spins.offset = linspace(-16e6,16e6,250);

% pi/2 pulse

Sim.direction = 'y';
Sim.angle = 90;
Sim.nframes = 30;

Spins.degrees = pulsespin(Spins,Sim);

% tau

Sim.length = tau;
Sim.nframes = 50;

Spins.degrees = evolution(Spins,Sim);

% first pi pulse

Sim.direction = 'y';
Sim.angle = 180;
Sim.nframes = 60;

Spins.degrees = pulsespin(Spins,Sim);

% time between pi pulse and pump pulse

Sim.length = tau+timing;
Sim.nframes = round(200*((timing+tau)/(T-timing)));

Spins.degrees = evolution(Spins,Sim);

% changes the precesion frequnecy for 10% of the spins according to the
% eecoupling

for i=200:250;
    Spins.offset(i) = Spins.offset(i)-offset;
end

% time between pump pulse and refocusing pulse

Sim.length = T-timing;
Sim.nframes = round(200*(T-timing)/(T+tau));

Spins.degrees = evolution(Spins,Sim);

% refocusing pulse

Sim.direction = 'y';
Sim.angle = 180;
Sim.nframes = 60;

Spins.degrees = pulsespin(Spins,Sim);

% time between refocussing pulse and the top of the refocussed echo

Sim.length = T;
Sim.nframes = 100;

Spins.degrees = evolution(Spins,Sim);


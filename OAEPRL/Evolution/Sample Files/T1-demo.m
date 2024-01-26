%% This simulation shows the effect of T1 following a pi pulse

clear all

%% Creates a matrix of spins along the z axis

for i=1:1
    Spins.degrees(i,:) = [0,0,1];
end

Spins.offset = linspace(0,0,1);
Spins.Tone = 0.0002;
Spins.Ttwo = 0.000002;

%% Sets the parameters and function for first pi pulse

Pulse.direction = 'y';
Pulse.angle = 180;

Spins.degrees = pulsespin(Spins,Pulse);

%% Sets the evolution parameters and function between pulses

Sim.length = 0.002;
Sim.nframes = 100;
Sim.delay = 0.1;

Spins.degrees = evolution(Spins,Sim);
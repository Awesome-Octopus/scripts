%% This simulation shows the effect of T2 following a pi/2 pulse

clear all

%% Creates a matrix of spins along the z axis

for i=1:5
    Spins.degrees(i,:) = [0,0,1];
end

Spins.offset = linspace(-50000,50000,5);
Spins.Tone = 0.005;
Spins.Ttwo = 0.0000115;

%% Sets the parameters and function for first pi/2 pulse

Pulse.direction = 'y';
Pulse.angle = 90;

Spins.degrees = pulsespin(Spins,Pulse);

%% Sets the evolution parameters and function between pulses

Sim.length = 0.00007;
Sim.nframes = 100;

Spins.degrees = evolution(Spins,Sim);
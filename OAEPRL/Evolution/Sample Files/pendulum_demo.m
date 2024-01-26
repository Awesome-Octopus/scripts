clear all

%% Creates a matrix of spins along the z axis

for i=1:15
    Spins.degrees(i,:) = [0,90,1];
end

%% Sets the evolution parameters and function between pulses

Spins.Tone = 10000;
Spins.Ttwo = 10000;
Spins.offset = linspace(51/60,65/60,15);

Sim.length = 180;
Sim.nframes = 1800

Spins.degrees = evolution(Spins,Sim);
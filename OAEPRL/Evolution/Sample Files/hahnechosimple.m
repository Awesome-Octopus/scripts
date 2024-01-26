%%  This imulation shows a simple Hahn Echo

clear all

%% Creates a matrix of spins along the z axis

for i=1:5
    Spins.degrees(i,:) = [0,0,1];
end

Spins.offset = linspace(-16e5,16e5,5);
Spins.Tone = 0.0002;
Spins.Ttwo = 0.000002;

%% Sets the parameters and function for first pi/2 pulse

Sim.direction = 'y';
Sim.angle = 90;
Sim.framerate = 20;
Sim.nframes = 30;

Spins.degrees = pulsespin(Spins,Sim);

%% Sets the evolution parameters and function between pulses

Sim.length = 0.0000002;
Sim.nframes = 50;

Spins.degrees = evolution(Spins,Sim);

%% Sets the parameters and function for pi pulse

Sim.direction = 'y';
Sim.angle = 180;
Sim.nframes = 60;

Spins.degrees = pulsespin(Spins,Sim);

%% Sets the evolution parameters and function after pulses

Sim.length = 0.0000004;
Sim.nframes = 100;

Spins.degrees = evolution(Spins,Sim);

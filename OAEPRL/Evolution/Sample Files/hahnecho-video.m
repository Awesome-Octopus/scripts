%%  This script shows the method for creating an animation from your simulation

clear all

%% Sets the parameters for the animation

Sim.video = 'avi'; % change this to 'gif' to get an animated GIF file
Sim.filename = '900500.avi'; % Change this to '900500.gif' if 'gif' is specified above
Sim.figuresize = [900 500];
Sim.numloops = 0;
Sim.framerate = 30;

%% Creates a matrix of spins along the z axis

for i=1:250
    Spins.degrees(i,:) = [0,0,1];
end

Spins.offset = linspace(-1600000,1600000,250);
Spins.Tone = 0.0002;
Spins.Ttwo = 0.00002;

%% Sets the parameters and function for first pi/2 pulse

Sim.direction = 'y';
Sim.angle = 90;
Sim.nframes = 45;

[Spins.degrees Sim] = pulsespin(Spins,Sim);

%% Sets the evolution parameters and function between pulses

Sim.length = 0.0000005;
Sim.nframes = 90;

[Spins.degrees Sim] = evolution(Spins,Sim);

%% Sets the parameters and function for first pi/2 pulse

Sim.direction = 'y';
Sim.angle = 180;
Sim.nframes = 90;

[Spins.degrees Sim] = pulsespin(Spins,Sim);

%% Sets the evolution parameters and function between pulses

Sim.length = 0.0000010;
Sim.nframes = 180;

[Spins.degrees Sim] = evolution(Spins,Sim);
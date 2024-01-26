clear all

%% Creates a matrix of spins along the z axis

Sim.video = 'y'

for i=1:5
    Spins.degrees(i,:) = [0,0,1];
end

%% Sets the parameters and function for first pi pulse

Sim.direction = 'y';
Sim.angle = 180;
Sim.delay = 0.002;

[Spins.degrees Sim.frame] = pulse(Spins,Sim);

%% Sets the evolution parameters and function between pulses

Spins.Tone = 0.0002;
Spins.Ttwo = 0.002;
Spins.offset = linspace(-160000,160000,5);

Sim.length = 0.00002;
Sim.inc = 0.000002;
Sim.delay = 0.1;

[Spins.degrees Sim.frame] = evolution(Spins,Sim);

%% Sets the parameters and function for first pi/2 pulse

Sim.direction = 'y';
Sim.angle = 90;
Sim.delay = 0.002;

[Spins.degrees Sim.frame] = pulse(Spins,Sim);

Sim.length = 0.0000005;
Sim.inc = 0.00000001;
Sim.delay = 0.05;

[Spins.degrees Sim.frame] = evolution(Spins,Sim);

%% Sets the parameters and function for pi pulse

Sim.direction = 'y';
Sim.angle = 180;
Sim.delay = 0.01;

[Spins.degrees Sim.frame] = pulse(Spins,Sim);

%% Sets the evolution parameters and function after pulses

Sim.length = 0.0000008;

[Spins.degrees Sim.frame] = evolution(Spins,Sim);
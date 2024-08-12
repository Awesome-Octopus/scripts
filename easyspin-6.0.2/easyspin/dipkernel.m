% dipkernel   Kernel matrix for DEER
%
%  K = dipkernel(t,r)
%  K = dipkernel(t,r,[g1 g2])
%
%  Provides the dipolar kernel matrix for DEER spectroscopy, assuming a
%  complete powder average and no orientation selection.
%
%  Inputs:
%    t ...       array of time values (microseconds)
%    r ...       array of distance values (nanometers)
%    [gA gB] ... the g values of the two spins, assumed 2.002319 if not
%                given
%
%  Output:
%    K ... kernel matrix where K(i,j) is K(t(i),r(j))
%
%  Example: Calculate the DEER signal from a distance distribution P(r)
%
%  r = linspace(0,6,1001);   % distance range
%  P = gaussian(r,3.5,0.3);  % distance distribution
%  t = linspace(-0.2,4,301); % time range
%  K = dipkernel(t,r);       % dipolar kernel
%  V = K*P(:);               % dipolar signal averaged over distribution
%  plot(t,V);
%  

% fastmotion  Computes linewidth parameters for fast-motion regime
%
%   lw = fastmotion(Sys,Field,tcorr)
%   lw = fastmotion(Sys,Field,tcorr,domain)
%   [lw,mI] = fastmotion(...)
%
%   Computes FWHM Lorentzian linewidths
%   for the fast motional regime according to
%
%       FWHM = A + B*mI + C*mI^2
%
%   (one nucleus) and a similar formula for systems
%   with more than one nucleus. If present, the contribution
%   of the nuclear quadrupole is also included.
%
%   System    spin system structure
%      g      principal values of g matrix
%      gFrame g Euler angles (optional) [radians]
%      Nucs   nuclear isotope(s), e.g. '14N' or '63Cu'
%      A      principal values of A matrix [MHz]
%      AFrame A Euler angles (optional) [radians]
%      Q      principal values of Q tensor [MHz]
%   Field     center magnetic field [mT]
%   tcorr     isotropic rotational correlation time [seconds]
%   domain    'feield' or 'freq', for linewidths
%
%   lw        all FWHM line widths [mT or MHz]
%   mI        mI quantum numbers for the lines, one line per row
%
%   Example:
%
%     System = struct('g',[2.0088 2.0064 2.0027],'Nucs','14N');
%     System.A = mt2mhz([7.59 5.95 31.76]/10);
%     [lw,mI] = fastmotion(System,350,1e-10)

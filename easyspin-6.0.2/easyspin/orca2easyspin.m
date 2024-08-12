% orca2easyspin   Import spin Hamiltonian parameters from ORCA
%
%  Sys = orca2easyspin(OrcaFileName)
%  Sys = orca2easyspin(OrcaFileName,HyperfineCutoff)
%
%  Loads the spin Hamiltonian parameters from the ORCA output file
%  given in OrcaFileName and returns them as an EasySpin spin system
%  structure Sys. If the output file contains multiple structures, Sys
%  is an array of spin system structures.
%
%  Input:
%    OrcaFileName     file name of the main ORCA output file
%    HyperfineCutoff  cutoff for hyperfine coupling (MHz)
%
%  Output:
%    Sys       spin system structure, or array of spin system structures
%    Sys.data  contains additional data read from the output file
%              (coordinates, charge, electric field gradients, etc)
%
%  Besides the main text-formatted output file, ORCA also generates an
%  additional file that contains atomic coordinates and calculated
%  properties such as g and A matrices, Q tensors, etc. This property file
%  is text-based and ends in _property.txt. Before ORCA 5, the property
%  file was binary and had extension .prop. orca2easyspin can read either
%  the main output file or the associated property file.
%
%  Examples:
%    Sys = orca2easyspin('nitroxide.out')   % all ORCA versions
%    Sys = orca2easyspin('nitroxide_property.txt')   % ORCA v5 and later
%    Sys = orca2easyspin('nitroxide.prop')   % ORCA prior to v5
%
%  If HyperfineCutoff (a single value, in MHz) is given, all nuclei with
%  hyperfine coupling equal or smaller than that value are omitted from
%  the spin system. If not given, it is set to zero, and all nuclei with
%  non-zero hyperfine coupling are included.
%
%  Example:
%    Sys = orca2easyspin('nitroxide.out',0.5)  % 0.5 MHz hyperfine cutoff

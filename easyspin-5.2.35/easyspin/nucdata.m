% nucdata  Nuclear spin data 
%
%   Spin = nucdata(Isotopes)
%   [Spin,gn] = nucdata(Isotopes)
%   [Spin,gn,qm] = nucdata(Isotopes)
%   [Spin,gn,qm,abund] = nucdata(Isotopes)
%
%   Returns nuclear spin, gyromagnetic
%   ratio, quadrupole moment and natural
%   abundance of one or several nuclei.
%
%   Isotopes is a string specifying the
%   nucleus, e.g. '1H', '13C', '63Cu', '191Ir'.
%   If Isotopes is a comma-separated list of nuclei
%   like '14N,14N,14N,1H,63Cu', vectors are returned.
%   Attention, case-sensitive!

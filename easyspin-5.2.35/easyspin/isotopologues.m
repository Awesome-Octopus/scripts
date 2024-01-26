% isotopologues   Generate list of isotopologues
%
%    isotopologues(NucList)
%    isotopologues(NucList,nEquiv)
%    isotopologues(NucList,nEquiv,Abundances)
%    isotopologues(NucList,nEquiv,Abundances,relThreshold)
%    isotopologues(Sys)
%    isotopologues(Sys,relThreshold)
%    out = ...
%
%    Computes all possible isotopologues of the given list of nuclei or
%    elements, including their abundances.
%
%    Input:
%      NucList       list of isotopes or elements
%                       e.g. 'Cu', '63Cu,N,N', 'H,H,C,C'
%      nEquiv        number of equivalent nuclei (default: 1 for each
%                       nucleus in NucList)
%                       e.g. [1 4] for 'Cu,N'
%      Abundances    cell array of nuclear abundances
%      relThreshold  isotopologue abundance threshold, relative to
%                       abundance of most abundant isotopologue
%                       (default 0.001)
%
%    out                 structure array containing a list of all isotopologues
%       out(k).Nucs      string with list of isotopes
%       out(k).Abund     overall absolute abundance
%       out(k).n         number of equivalent nuclei
%
%    If out is not requested, the list of isotopologues is displayed.

% ctafft  Cross-term averaged FFT 
%
%   FD = ctafft(TD,Averages)
%   FD = ctafft(TD,Averages,N)
%
%   Similar to a magnitude FFT, but removes
%   dead-time artefacts using cross-term
%   averging.
%   TD: Time-domain data. For matrices
%      TD, ctafft works along columns.
%   Averages: Starting indices for the FFTs
%     for TD(Averages(i):end). IF Averages is
%     a single number, it stands for 1:Averages.
%   N: Length of FFT. If omitted, N is set to
%     the length of TD (or its columns).

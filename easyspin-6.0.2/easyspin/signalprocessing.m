% signalprocessing    Signal translation and clean up for spidyan
%
%     x = signalprocessing(TimeAxis,RawSignal,FreqTranslation)
% 
%     TimeAxis   ... time axis of the input signal in microseconds
%     RawSignal   ... original signal, as returned by saffron and spidyan
%     FreqTranslation   ... vector with frequencies in GHz for translation
%
%     out:
%       x          ... signal(s) after frequency translation

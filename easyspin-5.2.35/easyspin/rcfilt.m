% rcfilt  Apply an RC filter to a signal
%
%  yFilt = rcfilt(y,SampTime,TimeConstant)
%  yFilt = rcfilt(y,SampTime,TimeConstant,'up')
%  yFilt = rcfilt(y,SampTime,TimeConstant,'down')
%
%  Filters a spectrum using a RC low-pass
%  filter as built into cw spectrometers
%  to remove high-frequency noise.
%
%  Input:
%  - y              unfiltered spectrum
%  - SampleTime     sampling time
%  - TimeConstant   time constant of the filter
%  - 'up' or 'down' defines the direction of the
%       field sweep. If omitted, 'up' is assumed.
%
%    SampleTime and TimeConstant must have the same
%    units (s, ms, us, ...)
%
%  Output:
%  - yFilt          filtered spectrum
%
%  For matrices, rcfilt operates along columns.

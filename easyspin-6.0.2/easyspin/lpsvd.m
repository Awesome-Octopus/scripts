% Linear Prediction Singular Value Decomposition
%
%  predictedSpectrum = lpsvd(Spectrum, Time, Method, Order)
% [predictedSpectrum, PredictionParameters] = lpsvd(...)
%
% Performs Linear Prediction SVD using a damped exponential model:
% y = amp * exp(1i*phase) * exp(time * (1i*2*pi*freq - damp) );
%
% Spectrum - The spectrum to be predicted
% Time - The corresponding time vector for the spectrum, necessary to
%        predict the proper phase of the signal.
%
% Method - the method string input determines the LPSVD algorithm used,
%          if not provided it will default to 'ss'
%
% Methods:
% 'kt' Based on:
% Kumaresan, R.;Tufts, D.W.; IEEE Trans. Acoust. Speech Signal ASSP-30 833 (1982)
%
% 'ss'  state-space method:
% Kung, S.Y.; Arun, K.S.; Bhaskar Rao, D.V.; J. Opt. Soc. Am. 73, 1799 (1983)
% Barkhuijsen, H.; De Beer, R.; Van Ormondt, D.; J. Mag. Reson. 73, 553 (1987)
%
% 'tls' Hankel Total least squares method:
% Van Huffel, S.;Chen, H.; Decanniere, C.; Van Hecke, P.; J. Mag. Reson.,
% Series A 110, 228 (1994)
%
%
% Order - The number of sinusoides to attempt to fit to the data if no order
%         is provided it will be automatically estimated.
%
% Model Order estimated as per:
% Wax, M.; Kailath, T.; IEEE Trans. Acoust. Speech Signal ASSP-39 387 (1985)
% 'mdl'  minimum description length
% 'aic' Akaike information protocol
% however these methods are known to underestimate the number of components

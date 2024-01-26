% apowin  Apodization windows 
%
%    w = apowin(Type,n)
%    w = apowin(Type,n,alpha)
%
%    Returns an apodization window. n is the number
%    of points. The string Type specifies the type
%    and can be
%
%      'bla'    Blackman
%      'bar'    Bartlett
%      'con'    Connes
%      'cos'    Cosine
%      'ham'    Hamming
%      'han'    Hann (also called Hanning)
%      'wel'    Welch
%
%   The following three windows need the parameter
%   alpha. Reasonable ranges for alpha are given.
%
%      'exp'    Exponential    2 to 6
%      'gau'    Gaussian       0.6 to 1.2
%      'kai'    Kaiser         3 to 9
%
%   A '+' ('-') appended to Type indicates that only the
%   right (left) half of the window should be constructed.
%
%      'ham'    symmetric (-1 <= x <= 1, n points)
%      'ham+'   right side only (0 <= x <= 1, n points)
%      'ham-'   left side only (-1 <= x <= 0, n points)

% basecorr  Baseline correction
%
%   datacorr = basecorr(data,dim,n)
%   datacorr = basecorr(data,dim,n,region)
%   [datacorr,baseline] = basecorr(...)
%
%   Performs polynomial baseline correction.
%
%   Inputs:
%     data       data to baseline correct, 1D or 2D
%     dim        dimension along which to baseline correct; if [], then
%                  a 2D function is fitted.
%     n          polynomial order, scalar for 1D data or row/column-wise 2D data,
%                and 2-element array for 2D; between 0 and 6
%     region     logic array, providing a mask for the region to include in
%                the baseline fit; true for each point to be included, and
%                false for each point to be excluded
%
%   Output:
%     datacorr   baseline-corrected data
%     baseline   fitted baseline
%
%   Examples:
%      z = basecorr(z,1,2);         % 2nd-order along dimension 1
%
%      region = x<10 | x>50;        % include x regions <10 and >50
%      z = basecorr(z,1,2,region);
%
%      z = basecorr(z,3,[]);        % two-dimensional 3rd-order surface
%

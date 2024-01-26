function APT=get_APT_kernel(ndip)
%
% get correct APT kernel
%

pot=0; rem=ndip; % find minimum power of 2 which is still larger than number of data points 
while rem>1,
    pot=pot+1;
    rem=rem/2;
end;
ksize=2^pot; % actual kernel size
p_string=sprintf('%d',ksize); % display
fname=sprintf('%s%s','kernel',p_string); % generate filename
load(fname);
APT.kernel=base;
APT.crosstalk=crosstalk;
APT.norm=tnorm;
APT.ny=ny;
APT.t=t;

function bckg=fit_bckg_std_alone(texp,t_fit,td_fit,back_model,hom_dim)
%
% Fits background decay according to the model selected by
% back_model, which defaults to 2
% 2  three-dimensional, exp(-k*t)
% 3  polynomial with order hom_dim (hom_dim defaults to 3)

if nargin<4,
    back_model=2;
end;

if nargin<5,
    hom_dim=3;
end;

        
poly=polyfit(t_fit,log(td_fit),1); % linear fit of logarithm
switch back_model,
    case 2,
		bckg=exp(polyval(poly,texp)); % background is exponential of that
    case 3,
        poly=polyfit(t_fit,log(td_fit),hom_dim);
        bckg=exp(polyval(poly,texp));
end;

bckg=real(bckg);

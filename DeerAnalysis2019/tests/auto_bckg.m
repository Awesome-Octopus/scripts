function [ff,bckg]=auto_bckg(texp,vexp)

ff=[];
bckg=[];

nexp=length(texp);
vexp=vexp/max(real(vexp));

APT=get_APT_kernel(nexp);
base=APT.kernel;
crosstalk=APT.crosstalk;
trcnorm=APT.norm;
t=APT.t;


% Adaptive background correction
nfa=round(0.1*nexp);
if nfa<1, nfa=1; end;
nfe=round(0.6*nexp);
if nfe<5, nfe=5; end;
merit=zeros(1,nfe-nfa);
%     aptmat=zeros(nfe-nfa,149);
best_merit=1e6;
for nofitp0=nfa:nfe,
%         disp(sprintf('%s%i','Testing background start at: ',nofitp0));
	t_fit=texp(nofitp0:length(texp)); % time window of baseline region
	td_fit=vexp(nofitp0:length(vexp)); % experimental data in this window
	
	% Background fit    
	td_poly=fit_bckg_std_alone(texp,t_fit,td_fit);
    
	td_exp2=vexp-td_poly; % subtract background
	td_exp2=td_exp2./td_poly; % divide by background, eqn [13]
	%dipevo=td_exp2/cfac2;
	dipevo=td_exp2/max(td_exp2); % normalize
    [m,n]=size(base); % size of kernel
	spc2=zeros(1,m); % initialize distribution
    td=zeros(1,n);
    td(1:length(texp))=dipevo;
	tdx=td.*t; % eqn [21]
	for k=1:m, % sum in eqn [21]
      spc2(k)=spc2(k)+sum(base(k,:).*tdx)/trcnorm(k); 
	end;
	spc3=crosstalk\spc2'; % crosstalk correction, eqn [22]
    merit(nofitp0-nfa+1)=sum(abs(spc3(1:3)));
    tact=texp(nofitp0);
    % fprintf(1,'%s%d%s%6.4f\n','Optimizing background fit range, Start: ',tact,' ns, Figure of merit: ',abs(spc3(1)));
    if merit<best_merit,
        ff=dipevo;
        bckg=td_poly;
    end;
end;

wdd0=52.23;

rax=1.5:0.01:8;
rmean=3;
sigr=0.01*sqrt(2);

arg=(rax-rmean)/sigr;
distr=exp(-arg.^2);
distr=distr/sum(distr);

t=0:0.008:2.048;
ff=zeros(size(t));

tic,
for k1=1:length(rax),
    wdd=wdd0/(rax(k1)^3);
    p=distr(k1);
    for kx=0:1000,
        x=kx/1000;
        w=2*pi*(3*x^2-1)*wdd;
        ff=ff+cos(w*t)*p;
    end;
end
toc,
ff=ff/max(ff);

sigr=0.15*sqrt(2);

arg=(rax-rmean)/sigr;
distr1=exp(-arg.^2);
distr1=distr1/sum(distr1);
ff1=zeros(size(t));

tic,
for k1=1:length(rax),
    wdd=wdd0/(rax(k1)^3);
    p=distr1(k1);
    for kx=0:1000,
        x=kx/1000;
        w=2*pi*(3*x^2-1)*wdd;
        ff1=ff1+cos(w*t)*p;
    end;
end
toc,
ff1=ff1/max(ff1);

mom=moment_analysis_vec(rax,distr1);
fprintf(1,'Gaussian mean : %6.4f nm, width: %6.4f nm\n',mom(1),real(sqrt(mom(2))));

ff2=zeros(size(t));

rmean2=rmean^2/2.8729;
sigr2=0.1571;
%rmean2=rmean;
%sigr2=sigr;
kappa=4;
distr2=zeros(size(rax));
for k=1:length(rax),
    if rax(k)<=rmean2,
        distr2(k)=exp(-(rmean2-rax(k))^2/(sigr2^2*kappa));
    else
        distr2(k)=exp(-kappa*(rax(k)-rmean2)^2/sigr2^2);
    end;
end;
distr2=distr2/sum(distr2);
tic,
for k1=1:length(rax),
    wdd=wdd0/(rax(k1)^3);
    p=distr2(k1);
    for kx=0:1000,
        x=kx/1000;
        w=2*pi*(3*x^2-1)*wdd;
        ff2=ff2+cos(w*t)*p;
    end;
end
toc,
ff2=ff2/max(ff2);
mom=moment_analysis_vec(rax,distr2);
fprintf(1,'Asymmetric mean : %6.4f nm, width: %6.4f nm\n',mom(1),real(sqrt(mom(2))));

figure(11); clf;
plot(t,ff,'k:');
hold on;
plot(t,ff1,'r--');
plot(t,ff2,'b');
set(gca,'FontSize',14);
axis([0,2,-0.25,1.05])

figure(12); clf;
% plot(rax,distr,'k');
hold on
plot(rax,distr1,'r--');
plot(rax,distr2,'b');
set(gca,'FontSize',14);
axis([2,4,-0.1*max(distr2),1.1*max(distr2)]);


function [min_rmsd,regpar,rmean,sigr,rmean_found,sigr_found,ff,ff1,best_sim,deer,bckg]=reanalyze_DEER(testpair,noiselev,tcut,tau,regpar_vec)

if nargin<5,
    regpar_vec=[0.001,0.01,0.1,1,10,100,1000,10000,100000];
end;
if nargin<4,
    tau=5;
end;
if nargin<3,
    tcut=9;
end;
if nargin<2,
    noiselev=0.005;
end;
if nargin<1,
    testpair=7;
end;

switch testpair
    case 1
        basname='[2WSX](A){1}179_[2WSX](A){1}378';
    case 2
        basname='[2WSX](A){1}318_[2WSX](A){1}249';
    case 3
        basname='[2WSX](A){1}33_[2WSX](A){1}70';
    case 4
        basname='[2WSX](A){1}344_[2WSX](A){1}374';
    case 5
        basname='[2WSX](A){1}374_[2WSX](A){1}433';
    case 6
        basname='[2WSX](A){1}433_[2WSX](A){1}50';
    case 7
        basname='[2WSX](A){1}472_[2WSX](A){1}249';
end;

clc;

dname=[basname '_distr.dat'];
data=load(dname);
rax=data(:,1);
distr0=data(:,2);

mom=moment_analysis_vec(rax,distr0);
rmean=mom(1);
sigr=real(sqrt(mom(2)));

fname=[basname '_fit.dat'];
data=load(fname);
t=data(:,1);
deer0=data(:,2);

if length(t)>999,
    t=t(1:999);
    deer0=deer0(1:999);
end;

if tcut<t(end),
    [mi,poi]=min(abs(t-tcut));
    t=t(1:poi);
    deer0=deer0(1:poi);
end;

ff=deer0-0.6;
ff=ff/max(ff);
dmp=exp(-t/tau);
ff_bckg=deer0.*dmp;
deer=ff_bckg+noiselev*randn(size(ff_bckg));
dat2=[1000*t,deer];
save('test.dat','dat2','-ascii');

[ff1,bckg]=auto_bckg(t',deer');

eta_vec=zeros(1,length(regpar_vec));
rho_vec=zeros(1,length(regpar_vec));
rmsd_vec=zeros(1,length(regpar_vec));
rmean_vec=zeros(1,length(regpar_vec));
sigr_vec=zeros(1,length(regpar_vec));

min_rmsd=1e6;
best_k=1;
for k=1:length(regpar_vec),
    [r,distr,sim,rho,eta]=get_Tikhonov(t',ff1,regpar_vec(k));
    [mi,poi]=min(abs(r-7));
    r=r(1:poi);
    distr=distr(1:poi);
    [mi,poi]=min(abs(r-1.75));
    r=r(poi:end);
    distr=distr(poi:end);
    distr1=interp1(rax,distr0,r,'pchip',0);
    sc1=sum(distr1.*distr1)/sum(distr1.*distr);
    sim=sim'-0.99;
    sc2=sum(ff.*ff)/sum(sim.*ff);
    rms_distr=sqrt(sum((distr1-sc1*distr).^2)/length(distr));
    rmsd_vec(k)=rms_distr;
    eta_vec(k)=eta;
    rho_vec(k)=rho;
    mom=moment_analysis_vec(r,distr);
    rmean_vec(k)=mom(1);
    sigr_vec(k)=real(sqrt(mom(2)));
    if rms_distr<min_rmsd,
        min_rmsd=rms_distr;
        best_k=k;
        best_sim=sc2*sim;
        best_distr=sc1*distr;
    end;
end;



figure(1);
clf;
plot(rax,distr0,'k');
hold on
plot(r,best_distr,'r');
set(gca,'FontSize',14);
axis([min(r),max(r),-0.1*max(distr0),1.1*max(distr0)]);

figure(2); clf;
plot(t,ff,'k');
hold on;
plot(t,ff1,'g');
plot(t,best_sim,'r');

figure(3); clf;
plot(t,deer,'k');
hold on;
plot(t,bckg,'r');
set(gca,'FontSize',14);
axis([0,max(t),min(deer)-0.05,1.05]);

figure(4); clf;
plot(rho_vec,eta_vec,'ko');
hold on;
plot(rho_vec(best_k),eta_vec(best_k),'ko','MarkerFaceColor','r');

figure(5); clf;
plot(log10(regpar_vec),rmean_vec,'ko');
hold on;
plot([log10(regpar_vec(1)),log10(regpar_vec(end))],[rmean,rmean],'g:');

figure(6); clf;
plot(log10(regpar_vec),sigr_vec,'ko');
hold on;
plot([log10(regpar_vec(1)),log10(regpar_vec(end))],[sigr,sigr],'g:');

figure(7); clf;
plot(log10(regpar_vec),rmsd_vec,'ko');
hold on;
plot(log10(regpar_vec(best_k)),rmsd_vec(best_k),'ko','MarkerFaceColor','r');

fprintf(1,'--- At a noise level of %6.3f ---\n',noiselev); 
fprintf(1,'Best r.m.s.d. agreement of distance distribution at regularization parameter: %g\n',regpar_vec(best_k)); 
fprintf(1,'Best r.m.s.d. : %9.6f\n',min_rmsd); 

mom=moment_analysis_vec(rax,distr0);
rmean=mom(1);
sigr=real(sqrt(mom(2)));
fprintf(1,'Actual mean distance    : %4.2f Å\n',mom(1)); 
fprintf(1,'Actual width            : %4.2f Å\n',real(sqrt(mom(2)))); 
mom=moment_analysis_vec(r,best_distr);
rmean_found=mom(1);
sigr_found=real(sqrt(mom(2)));
fprintf(1,'Determined mean distance: %4.2f Å\n',mom(1)); 
fprintf(1,'Determined width        : %4.2f Å\n',real(sqrt(mom(2)))); 

regpar=regpar_vec(best_k);


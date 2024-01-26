trials=250;
testpair=7;
tau=5;
regpar=10;
tau1=0.2;

tcut_vec=1:0.25:2.5;
T2=1.5;
dmp=exp(-2*(tcut_vec+2*tau1)/T2);
dmp=dmp./sqrt(tcut_vec); % same measurement time, more points => less acquistions, sqrt less S/N
dmp=dmp/dmp(1);
noiselev0=0.005;

rmean_vs_tcut=zeros(1,length(tcut_vec));
sigr_vs_tcut=zeros(1,length(tcut_vec));
rmsd_vs_tcut=zeros(1,length(tcut_vec));
rmean_sdev_vs_tcut=zeros(1,length(tcut_vec));
sigr_sdev_vs_tcut=zeros(1,length(tcut_vec));
rmsd_sdev_vs_tcut=zeros(1,length(tcut_vec));

tic,
for k=1:length(tcut_vec),
    tcut=tcut_vec(k);
    rmean_vec=zeros(1,trials);
    sigr_vec=zeros(1,trials);
    regpar_vec=zeros(1,trials);
    rmsd_vec=zeros(1,trials);
    noiselev=noiselev0/dmp(k);
    for kk=1:trials,
        [min_rmsd,regpar,rmean,sigr,rmean_found,sigr_found,ff,ff1,best_sim,deer,bckg]=reanalyze_DEER(testpair,noiselev,tcut,tau,regpar);
        rmean_vec(kk)=rmean_found;
        sigr_vec(kk)=sigr_found;
        rmsd_vec(kk)=min_rmsd;
    end;
    rmean_vs_tcut(k)=mean(rmean_vec);
    rmean_sdev_vs_tcut(k)=std(rmean_vec);
    sigr_vs_tcut(k)=mean(sigr_vec);
    sigr_sdev_vs_tcut(k)=std(sigr_vec);
    rmsd_vs_tcut(k)=sqrt(sum(rmsd_vec.^2)/length(rmsd_vec));
    rmsd_sdev_vs_tcut(k)=sqrt(std(rmsd_vec.^2));
end;

toc,

figure(11); clf;
errorbar(tcut_vec,rmean_vs_tcut,2*rmean_sdev_vs_tcut,'ko');
hold on
plot([min(tcut_vec),max(tcut_vec)],[rmean,rmean],'r:');

figure(12); clf;
errorbar(tcut_vec,sigr_vs_tcut,2*sigr_sdev_vs_tcut,'ko');
hold on;
plot([min(tcut_vec),max(tcut_vec)],[sigr,sigr],'r:');

figure(13); clf;
errorbar(tcut_vec,rmsd_vs_tcut,2*rmsd_sdev_vs_tcut,'ko');

save pair_2_scan_tcut noiselev0 dmp tcut_vec tau regpar rmean_vs_tcut rmean_sdev_vs_tcut sigr_vs_tcut sigr_sdev_vs_tcut rmsd_vs_tcut
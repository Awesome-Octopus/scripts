function [axis_vec,md] = plot_exp_dist(deer_basname,rax)

data = load(sprintf('%s_distr.dat',deer_basname));
rexp = 10*data(:,1);
distr_exp = data(:,2);
distr_lb = data(:,3);
distr_ub = data(:,4);
distr_exp = interp1(rexp,distr_exp,rax,'pchip',0);
scal = sum(distr_exp);
distr_exp = distr_exp/scal;
distr_lb = interp1(rexp,distr_lb,rax,'pchip',0);
distr_lb = distr_lb/scal;
distr_ub = interp1(rexp,distr_ub,rax,'pchip',0);
distr_ub = distr_ub/scal;
for k = 1:length(rax)
    plot([rax(k),rax(k)],[distr_lb(k),distr_ub(k)],'Color',[0.5,0.5,0.5]);
end
plot(rax,distr_exp,'k');
ma = max(distr_ub);
axis_vec = [0,rmax,-0.1*ma,1.05*ma];
axis(axis_vec);
xlabel('distance r [Å]');
ylabel('P(r)');

data = load(sprintf('%s_bckg.dat',deer_basname));
md = 1-data(1,3);
title(sprintf('m.d. %4.2f',md));
data = load(sprintf('%s_fit.dat',deer_basname));
tmax = data(end,1);
sc = (tmax/2)^(1/3);
r0 = 15;
r1 = 30*sc;
skip = false;
if r1 > rmax
    r1 = rmax;
    skip = true;
end
plot([r0,r1],[-0.04*ma,-0.04*ma],'Color',[0,0.6,0],'LineWidth',4);
if ~skip
    r2 = 40*sc;
    if r2 > rmax
        r2 = rmax;
        skip = true;
    end
    plot([r1,r2],[-0.04*ma,-0.04*ma],'Color',[0.8,0.8,0],'LineWidth',4);
end
if ~skip
    r3 = 50*sc;
    if r3 > rmax
        r3 = rmax;
        skip = true;
    end
    plot([r2,r3],[-0.04*ma,-0.04*ma],'Color',[0.8,0.6,0],'LineWidth',4);
end
if ~skip
    r4 = 60*sc;
    if r4 > rmax
        r4 = rmax;
    end
    plot([r3,r4],[-0.04*ma,-0.04*ma],'Color',[0.6,0,0],'LineWidth',4);
end

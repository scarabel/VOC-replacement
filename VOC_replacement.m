% Copyright (c) 2021 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE.txt for 
% licensing information. 
% 
% If using this code, please cite 
% Wu, Scarabel, Majeed, Bragazzi, Orbiski, The impact of public health
% interventions on delaying and mitigating against replacement by
% SARS-CoV-2 variants of concern, submitted (2021)
% 
%% Script VOC_replacement.m
% Generates the figures for the main text of the reference manuscript

clear
close all

saveplots = 0; % activate for saving figures
colorcode = lines(10);

% Initial parameter definition
y0 = 1000; % initial data
mean_gt = 7; % mean generation time
gamma = 1/mean_gt; % removal rate

Init_ratio = 0.01; % initial ratio between variant and resident
rel_transm = 1.5; % relative transmissibility of variant vs resident
R0 = 0.8; % control reproduction number of resident

r = gamma*(R0-1); % growth rate
Td = log(2)/r; % doubling time

% Figure 1: Plot of total cases
time=0:200;
cases_no_VOC = y0*exp(r*time);
total_cases = Init_ratio*exp(gamma*R0*(rel_transm-1)*time)+1.*cases_no_VOC/(1+Init_ratio);

figure(10); hold on
set(gca,'defaultAxesColorOrder',[colorcode(1,:); colorcode(2,:)]);
yyaxis left
plot(time,total_cases,'Linewidth',2,'Color',colorcode(1,:),'HandleVisibility', 'on')
plot(time,cases_no_VOC,'Linewidth',2,'Color',colorcode(1,:),'LineStyle','--','HandleVisibility', 'on')
legend('with variant','without variant','Location','NorthWest')
%title('Expected reported cases')
xlabel('Days')
ylabel('Daily reported cases')

% Figure 1: Plot of frequency
ratio = Init_ratio*exp(gamma*R0*(rel_transm-1)*time);

yyaxis right
plot(time,1./(1+1./ratio),'Linewidth', 2,'Color',colorcode(2,:),'HandleVisibility', 'off')
ylabel('Variant frequency')
set(gca,'FontSize',14)
set(gca,'YGrid','on','XGrid','on')

if saveplots
    set(gcf, 'Position', get(0, 'Screensize'));
    savefig(['exp-increase-R08']);
    saveas(gcf,['exp-increase-R08.png']);
end

%% Figure 2: Time to replacement as a function of the generation time

gen_time_vec = 1:0.2:10;
gamma_vec = 1./gen_time_vec';
R0_vec = [0.8;1.1];

Time_to_surpass = zeros(length(gamma_vec),length(R0_vec));
legendlist={};
figure(2); clf
subplot(1,2,1); hold on
for ind_R0=1:length(R0_vec)
    Time_to_surpass(:,ind_R0) = 1./(gamma_vec*R0_vec(ind_R0)*(rel_transm-1))*log(1/Init_ratio);
    plot(gen_time_vec,Time_to_surpass(:,ind_R0),'LineWidth',2,'Color',colorcode(ind_R0,:));
    legendlist{end+1} = ['$R_w=$',num2str(R0_vec(ind_R0))];
end
xlabel('Mean generation time (days)')
ylabel('Time to replacement (days)')
legend(legendlist,'Location','NorthWest','Interpreter','latex')
set(gca,'FontSize',14)
set(gca,'YGrid','on','XGrid','on')

%% Figure 2: Epidemic doubling time depending on the initial ratio

R0_vec = [0.8;1.1];
l0_vec = [0:0.01:1]';

Doubl_time = [0*l0_vec,0*l0_vec];
Doubl_time(1,:) = log(2)./(gamma*(R0_vec'*rel_transm-1));

Rlow=0.8;
lcrit = (1-Rlow)/(Rlow*rel_transm-1);
t0=1/(gamma*R0_vec(ind_R0)*(rel_transm-1))*log(lcrit/l0_vec(2));

subplot(1,2,2); hold on
legendlist = {};
for ind_R0 =1:length(R0_vec)
    for jj = 2:length(l0_vec)
        ysol = fsolve(@(x) l0_vec(jj)*exp(gamma*(R0_vec(ind_R0)*rel_transm-1)*x) + exp(gamma*(R0_vec(ind_R0)-1)*x) - 2*(1+l0_vec(jj)), 10*Doubl_time(jj-1,ind_R0));
        Doubl_time(jj,ind_R0) = ysol;
    end
    
    if R0_vec(ind_R0)<1
        Doubl_time(1,ind_R0)=NaN;
    else
        Doubl_time(1,ind_R0)=log(2)/(gamma*(R0_vec(ind_R0)-1));
    end
    
    plot(l0_vec,Doubl_time(:,ind_R0),'LineWidth',2,'Color',colorcode(ind_R0,:));
    legendlist{end+1} = ['$R_w=$',num2str(R0_vec(ind_R0))];
end

xlabel('Initial ratio of new variant over resident strain')
ylabel('Total doubling time (days)')
%title(['Mean generation time ',num2str(mean_gt),' days']);
legend(legendlist,'Location','NorthEast','Interpreter','latex')
axis([0 1 0 100])
set(gca,'FontSize',14)
set(gca,'YGrid','on','XGrid','on')

if saveplots
    set(gcf, 'Position', get(0, 'Screensize'));
    savefig(['time-to-replacement-double']); 
    saveas(gcf,['time-to-replacement-double.png']); 
end

%% Figure 3: Epidemic doubling time vs doubling time of resident strain

R0_vec = [0.7:0.01:1.5]'; 
i1 = find(R0_vec<1,1,'last');

% Varying generation time
l0_fixed = 0.01;
mean_gt_vec = [5:7]';
gamma_vec = 1./mean_gt_vec;
Doubl_time = zeros(length(R0_vec),length(gamma_vec));
Doubl_time(end,:) = log(2)./(gamma_vec*(R0_vec(end)-1));

figure(3); clf
subplot(1,2,1); hold on
for jj=1:length(gamma_vec)
    gamma = gamma_vec(jj);
    for ind_R0 = length(R0_vec)-1:-1:1
        ysol = fsolve(@(x) l0_fixed*exp(gamma*(R0_vec(ind_R0)*rel_transm-1)*x) + exp(gamma*(R0_vec(ind_R0)-1)*x) - 2*(1+l0_fixed), Doubl_time(ind_R0+1,jj));
        Doubl_time(ind_R0,jj) = ysol;
    end

    plot(log(2)./(gamma*(R0_vec(1:i1)-1)),Doubl_time(1:i1,jj),'LineWidth',2,'Color',colorcode(jj,:)); hold on
    plot(log(2)./(gamma*(R0_vec(i1+1:end)-1)),Doubl_time(i1+1:end,jj),'LineWidth',2,'Color',colorcode(jj,:),'HandleVisibility','off'); hold on

end

xlabel('Doubling time of resident strain (days)')
ylabel('Total doubling time (days)')
%title(['Mean generation time ',num2str(mean_gt),' days, initial ratio l0 =',num2str(l0_fixed)]);
legend('5 days','6 days','7 days','Location','NorthWest')
axis([-100 100 0 200])
set(gca,'FontSize',14)

% Varying initial frequency
l0_vec = [0.1; 0.05; 0.03]; % 0.1
mean_gt = 7;
gamma = 1./mean_gt;
Doubl_time = zeros(length(R0_vec),length(l0_vec));
Doubl_time(end,:) = log(2)./(gamma*(R0_vec(end)-1));

% figure; hold on
subplot(1,2,2); hold on

for jj=1:length(l0_vec)
    l0_fixed = l0_vec(jj);
    for ind_R0 = length(R0_vec)-1:-1:1
        ysol = fsolve(@(x) l0_fixed*exp(gamma*(R0_vec(ind_R0)*rel_transm-1)*x) + exp(gamma*(R0_vec(ind_R0)-1)*x) - 2*(1+l0_fixed), Doubl_time(ind_R0+1,jj));
        Doubl_time(ind_R0,jj) = ysol;
    end

    plot(log(2)./(gamma*(R0_vec(1:i1)-1)),Doubl_time(1:i1,jj),'LineWidth',2,'Color',colorcode(jj,:)); hold on
    plot(log(2)./(gamma*(R0_vec(i1+1:end)-1)),Doubl_time(i1+1:end,jj),'LineWidth',2,'Color',colorcode(jj,:),'HandleVisibility','off'); hold on

end

xlabel('Doubling time of resident strain (days)')
ylabel('Total doubling time (days)')
%title(['Mean generation time ',num2str(mean_gt),' days, initial ratio l0 =',num2str(l0_fixed)]);
legend('10%','5%','3%','Interpreter','latex','Location','NorthWest')
axis([-100 100 0 200])
set(gca,'FontSize',14)

if saveplots
    set(gcf, 'Position', get(0, 'Screensize'));
    savefig('doubling-vs-doubling');
    saveas(gcf,'doubling-vs-doubling.png');
end


%% Figure 4-5: Scenario projections of reported cases and relative frequency

Delay = 7;
startDate = 0; 
ProjDate2 = startDate+15+Delay; 
ProjDate3 = startDate+30+Delay; 
ProjDateEnd = startDate+75; 

addci = 1; % to add interval projections

y0 = 1000;
rel_freq0 = y0/(1+1/Init_ratio);

alpha = 1.5; % relative transmissibility
alpha_l = 1.4;
alpha_u = 1.6;

Rlow = 0.8; % reproduction numbers
Rhigh = 1.1;
rlow = gamma*(Rlow-1);
rhigh = gamma*(Rhigh-1);

xproj1 = [startDate:ProjDateEnd]';
xproj2 = [ProjDate2:ProjDateEnd]';
xproj3 = [ProjDate3:ProjDateEnd]';

ym1 = y0*exp(rlow*(xproj1-xproj1(1)));

ym1_VOC = (Init_ratio*exp(gamma*Rlow*(alpha-1)*(xproj1-xproj1(1)))+1).*ym1/(1+Init_ratio);
yl1_VOC = (Init_ratio*exp(gamma*Rlow*(alpha_l-1)*(xproj1-xproj1(1)))+1).*ym1/(1+Init_ratio);
yu1_VOC = (Init_ratio*exp(gamma*Rlow*(alpha_u-1)*(xproj1-xproj1(1)))+1).*ym1/(1+Init_ratio);

% lower branch
figure(20); hold on
plot(xproj1,ym1,'Linewidth', 2,'Color',colorcode(1,:),'HandleVisibility', 'on')

figure(21); hold on
plot(xproj1,ym1_VOC,'Linewidth', 2,'Color',colorcode(1,:),'HandleVisibility', 'on')
if addci
    fill([xproj1;flip(xproj1)],[yl1_VOC;flip(yu1_VOC)],colorcode(1,:),'FaceAlpha',0.2,'HandleVisibility','off')
end

% medium branch
y0_2 = y0*exp(rlow*(ProjDate3-xproj1(1)));
y0_2l = y0_2;
y0_2u = y0_2;

y0_VOC_2 = (Init_ratio*exp(gamma*Rlow*(alpha-1)*(ProjDate3-xproj1(1)))+1).*y0_2/(1+Init_ratio);
y0_VOC_2l = (Init_ratio*exp(gamma*Rlow*(alpha_l-1)*(ProjDate3-xproj1(1)))+1).*y0_2l/(1+Init_ratio);
y0_VOC_2u = (Init_ratio*exp(gamma*Rlow*(alpha_u-1)*(ProjDate3-xproj1(1)))+1).*y0_2u/(1+Init_ratio);

Freq0_2 = Init_ratio*exp((gamma+rlow)*(alpha-1)*(ProjDate3-xproj1(1)));
Freq0_2l = Init_ratio*exp((gamma+rlow)*(alpha_l-1)*(ProjDate3-xproj1(1)));
Freq0_2u = Init_ratio*exp((gamma+rlow)*(alpha_u-1)*(ProjDate3-xproj1(1)));

figure(20)
ym2 = y0_2*exp(rhigh*(xproj3-xproj3(1)));
yl2 = ym2; yu2 = ym2;

ym2_VOC = (Freq0_2*exp(gamma*Rhigh*(alpha-1)*(xproj3-xproj3(1)))+1).*ym2/(1+Freq0_2)*y0_VOC_2/y0_2;
yl2_VOC = (Freq0_2l*exp(gamma*Rhigh*(alpha_l-1)*(xproj3-xproj3(1)))+1).*yl2/(1+Freq0_2l)*y0_VOC_2l/y0_2l;
yu2_VOC = (Freq0_2u*exp(gamma*Rhigh*(alpha_u-1)*(xproj3-xproj3(1)))+1).*yu2/(1+Freq0_2u)*y0_VOC_2u/y0_2u;

plot(xproj3,ym2,'Linewidth', 2,'Color',colorcode(3,:),'HandleVisibility','on')

figure(21)
plot(xproj3,ym2_VOC,'Linewidth', 2,'Color',colorcode(3,:),'HandleVisibility', 'on')
if addci
    fill([xproj3;flip(xproj3)],[yl2_VOC;...
        flip(yu2_VOC)],colorcode(3,:),'FaceAlpha',0.2,'HandleVisibility','off')
end

% upper branch
y0_3 = y0*exp(rlow*(ProjDate2-xproj1(1)));
y0_3l = y0_3;
y0_3u = y0_3;

y0_VOC_3 = (Init_ratio*exp(gamma*Rlow*(alpha-1)*(ProjDate2-xproj1(1)))+1).*y0_3/(1+Init_ratio);
y0_VOC_3l = (Init_ratio*exp(gamma*Rlow*(alpha_l-1)*(ProjDate2-xproj1(1)))+1).*y0_3l/(1+Init_ratio);
y0_VOC_3u = (Init_ratio*exp(gamma*Rlow*(alpha_u-1)*(ProjDate2-xproj1(1)))+1).*y0_3u/(1+Init_ratio);

Freq0_3 = Init_ratio*exp((gamma+rlow)*(alpha-1)*(ProjDate2-xproj1(1)));
Freq0_3l = Init_ratio*exp((gamma+rlow)*(alpha_l-1)*(ProjDate2-xproj1(1)));
Freq0_3u = Init_ratio*exp((gamma+rlow)*(alpha_u-1)*(ProjDate2-xproj1(1)));

figure(20)
ym3 = y0_3*exp(rhigh*(xproj2-xproj2(1)));
yl3 = ym3; yu3 = ym3;

ym3_VOC = (Freq0_3*exp(gamma*Rhigh*(alpha-1)*(xproj2-xproj2(1)))+1).*ym3/(1+Freq0_3)*y0_VOC_3/y0_3;
yl3_VOC = (Freq0_3l*exp(gamma*Rhigh*(alpha_l-1)*(xproj2-xproj2(1)))+1).*yl3/(1+Freq0_3l)*y0_VOC_3l/y0_3l;
yu3_VOC = (Freq0_3u*exp(gamma*Rhigh*(alpha_u-1)*(xproj2-xproj2(1)))+1).*yu3/(1+Freq0_3u)*y0_VOC_3u/y0_3u;

plot(xproj2,ym3,'Linewidth', 2,'Color',colorcode(2,:),'HandleVisibility','on')

figure(21)
plot(xproj2,ym3_VOC,'Linewidth', 2,'Color',colorcode(2,:),'HandleVisibility', 'on')
if addci
    fill([xproj2;flip(xproj2)],[yl3_VOC;...
        flip(yu3_VOC)],colorcode(2,:),'FaceAlpha',0.2,'HandleVisibility','off')
end

%% Figure 4-5: Complete plot

figure(20)
limx = [startDate max(xproj2)]; 
set(gca,'Xlim',limx) 
ax20=gca;
ax20.YAxisLocation = 'right';

set(gca,'Ylim',[0,5000]);
limy = get(gca,'Ylim');
axis([0 ProjDateEnd 0 1500])
set(gca,'YGrid','on')
set(gca,'FontSize',14)

% plot intervention dates
plot([startDate-Delay,startDate-Delay], [0,limy(end)], 'k')
plot([startDate,startDate], [0,limy(end)], 'k:') % delay 
plot([ProjDate2-Delay,ProjDate2-Delay], [0,limy(end)], 'k','HandleVisibility', 'off') % relax interventions
plot([ProjDate2,ProjDate2], [0,limy(end)], 'k:','HandleVisibility', 'off') % delay
plot([ProjDate3-Delay,ProjDate3-Delay], [0,limy(end)], 'k','HandleVisibility', 'off') % relax interventions
plot([ProjDate3,ProjDate3], [0,limy(end)], 'k:','HandleVisibility', 'off') % delay

xlabel('Days')
ylabel('Confirmed cases')
legend('low projection','medium projection','high projection','intervention',[num2str(Delay),'-day delayed effect'],'Location','NorthWest')

if saveplots
    set(gcf, 'Position', get(0, 'Screensize'));
    savefig('proj_no_VOC')
    saveas(gca,'proj_no_VOC.png')
end

figure(21)
limx = [startDate max(xproj2)]; 
set(gca,'Xlim',limx) 
set(gca,'Ylim',limy);
axis([0 ProjDateEnd 0 1500])
set(gca,'YGrid','on')
set(gca,'FontSize',14)
ax21=gca;
ax21.YAxisLocation = 'right';

plot([startDate-Delay,startDate-Delay], [0,limy(end)], 'k')
plot([startDate,startDate], [0,limy(end)], 'k:') % delay 
plot([ProjDate2-Delay,ProjDate2-Delay], [0,limy(end)], 'k','HandleVisibility', 'off') % relax interventions
plot([ProjDate2,ProjDate2], [0,limy(end)], 'k:','HandleVisibility', 'off') % delay
plot([ProjDate3-Delay,ProjDate3-Delay], [0,limy(end)], 'k','HandleVisibility', 'off') % relax interventions
plot([ProjDate3,ProjDate3], [0,limy(end)], 'k:','HandleVisibility', 'off') % delay

xlabel('Days')
ylabel('Confirmed cases')
legend('low projection','medium projection','high projection','intervention',[num2str(Delay),'-day delayed effect'],'Location','NorthWest')

if saveplots
    set(gcf, 'Position', get(0, 'Screensize'));
    savefig('proj_VOC')
    saveas(gca,'proj_VOC.png')
end

%% Plot of relative frequency

figure(7); clf
hold on

% lower branch
Rel_freq_2 = Init_ratio*exp((gamma+rlow)*(alpha-1)*(xproj1-xproj1(1))); 
Rel_freq_2l = Init_ratio*exp((gamma+rlow)*(alpha_l-1)*(xproj1-xproj1(1)));
Rel_freq_2u = Init_ratio*exp((gamma+rlow)*(alpha_u-1)*(xproj1-xproj1(1)));
plot(xproj1,1./(1+1./Rel_freq_2),'Linewidth', 2,'Color',colorcode(1,:),'HandleVisibility', 'on')
% if addci
%     fill([xproj1;flip(xproj1)],[1./(1+1./Rel_freq_2l);...
%         flip(1./(1+1./Rel_freq_2u))],colorcode(1,:),'FaceAlpha',0.2,'HandleVisibility','off')
% end

% medium branch
Rel_freq_2 = Freq0_2*exp((gamma+rhigh)*(alpha-1)*(xproj3-xproj3(1)));
Rel_freq_2l = Freq0_2l*exp((gamma+rhigh)*(alpha_l-1)*(xproj3-xproj3(1)));
Rel_freq_2u = Freq0_2u*exp((gamma+rhigh)*(alpha_u-1)*(xproj3-xproj3(1)));

plot(xproj3,1./(1+1./Rel_freq_2),'Linewidth', 2,'Color',colorcode(3,:),'HandleVisibility','on')
% if addci
%     fill([xproj3;fliplr(xproj3)],[1./(1+1./Rel_freq_2l);...
%         fliplr(1./(1+1./Rel_freq_2u))],colorcode(3,:),'FaceAlpha',0.2,'HandleVisibility','off')
% end

% upper branch
Rel_freq_3 = Freq0_3*exp((gamma+rhigh)*(alpha-1)*(xproj2-xproj2(1)));
Rel_freq_3l = Freq0_3l*exp((gamma+rhigh)*(alpha_l-1)*(xproj2-xproj2(1)));
Rel_freq_3u = Freq0_3u*exp((gamma+rhigh)*(alpha_u-1)*(xproj2-xproj2(1)));

plot(xproj2,1./(1+1./Rel_freq_3),'Linewidth', 2,'Color',colorcode(2,:),'HandleVisibility','on')
% if addci
%     fill([xproj2;flip(xproj2)],[1./(1+1./Rel_freq_3l);...
%         flip(1./(1+1./Rel_freq_3u))],colorcode(2,:),'FaceAlpha',0.2,'HandleVisibility','off')
% end

limx = [startDate max(xproj2)]; 
set(gca,'Xlim',limx) 

limy = get(gca,'Ylim');
set(gca,'YGrid','on')
set(gca,'FontSize',14)
ax7=gca;
ax7.YAxisLocation = 'right';

plot([startDate,startDate-Delay], [0,limy(end)], 'k','HandleVisibility', 'on')
plot([startDate,startDate], [0,limy(end)], 'k:','HandleVisibility', 'on') % delay 
plot([ProjDate2-Delay,ProjDate2-Delay], [0,limy(end)], 'k','HandleVisibility', 'off') % relax interventions
plot([ProjDate2,ProjDate2], [0,limy(end)], 'k:','HandleVisibility', 'off') % delay
plot([ProjDate3-Delay,ProjDate3-Delay], [0,limy(end)], 'k','HandleVisibility', 'off') % relax interventions
plot([ProjDate3,ProjDate3], [0,limy(end)], 'k:','HandleVisibility', 'off') % delay

xlabel('Days')
ylabel('Relative frequency')
legend('low projection','medium projection','high projection','intervention',[num2str(Delay),'-day delayed effect'],'Location','NorthWest')

if saveplots
    set(gcf, 'Position', get(0, 'Screensize'));
    savefig('rel_freq')
    saveas(gca,'rel_freq.png')
end



f = figure();

% Han's
ICMS = [0.1275,0.1240,0.1247,0.1330,0.1541,0.1619];
ICMS_std_err = [0.0044,0.005,0.0061,0.0082,0.0034,0.0044];
bump = [0.1714,0.1714,0.1554,0.1554,0.1589,0.1589];
bump_std_err = [0.0034,0.0034,0.0021,0.0021,0.0017,0.0017];

f.Name = 'Han_16elecs_minRT_scatterSummary';
errorbar(bump,ICMS,bump_std_err,'horizontal','.','markersize',20,'color','k')
hold on
errorbar(bump,ICMS,ICMS_std_err,'vertical','.','markersize',20,'color','k')

hold on
plot([0,0.3],[0,0.3],'r--','linewidth',1.5)
xlabel('Bump RT (s)')
ylabel('ICMS RT(s)')
formatForLee(gcf)
set(gca,'fontsize',16)
xlim([0.1,0.2])
ylim([0.1,0.2])
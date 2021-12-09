 
 
BP1 = [0.492709,0.502044,0.485610,0.495895,0.559926,0.771988];
BP2 = [0.553139,0.543379,0.538420,0.551360,0.588097,0.782590];
BP3 = [0.685590,0.691129,0.701592,0.726736,0.770895,0.855659];
MF1 = [0.567877,0.571919,0.584604,0.603415,0.655931,0.801847];
MF2 = [0.501830,0.507873,0.519467,0.543055,0.587976,0.764079];
MF3 = [0.601295,0.597687,0.612666,0.634114,0.690806,0.800921];

xx = [0,0.2,0.4,0.6,0.8,1];
figure(1)
plot(xx,MF1,'b-o', 'LineWidth', 4, 'MarkerSize',12);
hold on
plot(xx, MF2,'r-x', 'LineWidth', 4, 'MarkerSize',14);
hold on
plot(xx, MF3,'g-*', 'LineWidth', 4, 'MarkerSize',14);
hold off

xticks([0,0.2,0.4,0.6,0.8,1])
xticklabels({'0', '20%', '40%', '60%', '80%', '100%'});
xlabel('Percent of edges')
ylabel('function prediction accuracy')
legend('BP 11-30', 'BP 31-100', 'BP 101-300', 'location', 'northwest');
set(gca, 'fontsize', 20)

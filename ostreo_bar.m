% bar charts with error bars for ostreo abundances along each transect
% NEED TO GET OI ON THERE TOO!!!


load ostreo_std.csv;
trans=[1;1;1;1;1;1;1;1;2;2;2;2;2;4;4;4;4;4;4;5;5;5;5;5;5;5];
A=find(trans==1);
B=find(trans==2);
D=find(trans==4);
E=find(trans==5);

figure(1)
bar(ostreo_std(A,1),ostreo_std(A,4),'w');hold on;
errorbar(ostreo_std(A,1),ostreo_std(A,4),ostreo_std(A,5),'.r');
% add figure labels and title
% saveas(gca, 'ostreo_Abar.eps','epsc');

figure(2)
bar(ostreo_std(B,1),ostreo_std(B,4),'w');hold on;
errorbar(ostreo_std(B,1),ostreo_std(B,4),ostreo_std(B,5),'.r');
% add figure labels and title
% saveas(gca, 'ostreo_Bbar.eps','epsc');

figure(3)
bar(ostreo_std(D,1),ostreo_std(D,4),'w');hold on;
errorbar(ostreo_std(D,1),ostreo_std(D,4),ostreo_std(D,5),'.r');
% add figure labels and title
% saveas(gca, 'ostreo_Dbar.eps','epsc');

figure(4)
bar(ostreo_std(E,1),ostreo_std(E,4),'w');hold on;
errorbar(ostreo_std(E,1),ostreo_std(E,4),ostreo_std(E,5),'.r');
% add figure labels and title
% saveas(gca, 'ostreo_Ebar.eps','epsc');


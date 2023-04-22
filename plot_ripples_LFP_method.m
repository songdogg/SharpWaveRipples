%% Script to generate LFP analysis figure

%% Run ripple_analysis.m until generation of binImage (break at line 82) to plot the following graphs

[pow, TFphase]  = morletTF(dsdata', (1:length(dsdata))./in.downFs, in.morletFreq);
basepow         = log(mean(pow(:,195500:196700),2));
pow             = log(pow) - repmat(basepow,1,size(pow,2));

plt             = 195000:200000;

subplot(5,1,1)
plot(plt,dsdata(plt))
ylabel({'Raw LFP','(µV)'},'FontWeight','bold'), xticklabels({})

subplot(5,1,2)
plot(plt,filterdata(plt))
ylabel({'Filtered LFP','(µV)'},'FontWeight','bold'), xticklabels({})

subplot(5,1,3)
plot(plt,zsc(plt)), hold on
plot(plt,zeros(length(plt),1),'m--','LineWidth',0.7)
plot(plt,zeros(length(plt),1)+3,'--','Color','#007200','LineWidth',0.7)
plot(zeros(length(-1:8),1)+198441,(-1:8)-0.5,'r:','LineWidth',1), 
fill([197441 197441 199441 199441],[-1.7 7.5 7.5 -1.7],'r','FaceAlpha',0.2,'EdgeColor','none'), hold off
text(198460,5.8,'IED \pm1s','Color','r','FontWeight','bold')
text(195050,5.8,'Minimum threshold','Color','m','FontWeight','bold'), text(195750,5.8,'Peak threshold','Color','#007200','FontWeight','bold')
ylabel({'LFP Amplitude','(Z-score)',''},'FontWeight','bold'), xticklabels({})

subplot(5,1,4)
imagesc(binImage(plt)'), colormap(gca,'gray')
yticks([]), ylabel({'LFP Amplitude','(Binary Image)',''},'FontWeight','bold')
xticks(0:500:5000), xticklabels({}), xlim([0 5000])

subplot(5,1,5)
contourf(1:length(plt), in.morletFreq, pow(:,plt), 50,'linestyle','none');
xlabel('Time (s)','FontWeight','bold'), xticklabels(xticks/in.downFs), ylabel({'Ripple Frequency','(Hz)'},'FontWeight','bold')
c = colorbar; c.Label.String = {'Power (dB)'}; c.Label.Color = 'w'; c.Label.FontSize = 9; c.Label.FontWeight = 'bold';
c.Location = 'east'; c.Color = 'w'; clim([0 3]); clear c

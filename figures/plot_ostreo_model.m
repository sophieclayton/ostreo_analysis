% need to modify this to plot all of the model runs on to a single figure
% and output at the right resolution/size for journal.
load data/ostreo.mat

ftype = 4;

if ftype == 1,
    hold on
figure;
    % for t=length(tstep);
    subplot(3,1,1);plot(y,OI(:,t),'k',S(kuro),O(kuro,1),'kx','LineWidth',1);axis([33.2 34.4 0 n000]);title('OI','FontSize',n)
    ylabel('copies ml^{-1}','FontSize',n);set(gca,'FontSize',n)
    legend('model','data','North');%legend BOXOFF
    %     hold on;plot(y,-(100/1.5)*y+2266.7,'r');hold off
    subplot(3,1,2);plot(y,OII(:,t),'k',S(kuro),O(kuro,2),'kx','LineWidth',1);axis([33.2 34.4 0 n000]);title('OII','FontSize',n)
    ylabel('copies ml^{-1}','FontSize',n);set(gca,'FontSize',n)
    %     hold on;plot(y,(100/1.5).*y-2200,'r');hold off
    subplot(3,1,3);plot(y,OII(:,t)./(OII(:,t)+OI(:,t)),'k',S(kuro),ostreo(kuro)./100,'kx','LineWidth',1);axis([33.2 34.4 0 1]);ylabel('clade ratio OII/OI','FontSize',n);
    xlabel('salinity','FontSize',n);set(gca,'FontSize',n)
    hold off
    
elseif ftype == 2,
    
figure;
    % for t=length(tstep);
    subplot(3,1,1);plot(y,OI_unetpos(:,t),'k',y,OI_unet0(:,t),'--k',y,OI_unetneg(:,t),':k',S(kuro),O(kuro,1),'kx','LineWidth',2,'MarkerSize',12);axis([33.2 34.4 0 18000]);title('OI','FontSize',n,'FontName','Helvetica')
    ylabel('Copies ml^{-1}','FontSize',n,'FontName','Helvetica');set(gca,'FontSize',n,'FontName','Helvetica','XTick',33.2:0.2:34.4)
    legend('\mu_{NET}>0','\mu_{NET}=0','\mu_{NET}<0','Kuroshio','North');%legend BOXOFF
    %     hold on;plot(y,-(100/1.5)*y+2266.7,'r');hold off
    subplot(3,1,2);plot(y,OII_unetpos(:,t),'k',y,OII_unet0(:,t),'--k',y,OII_unetneg(:,t),':k',S(kuro),O(kuro,2),'kx','LineWidth',2,'MarkerSize',12);axis([33.2 34.4 0 18000]);title('OII','FontSize',n,'FontName','Helvetica')
    ylabel('Copies ml^{-1}','FontSize',n,'FontName','Helvetica');set(gca,'FontSize',n,'FontName','Helvetica','XTick',33.2:0.2:34.4)
    %     hold on;plot(y,(100/1.5).*y-2200,'r');hold off
    subplot(3,1,3);plot(y,OII_unetpos(:,t)./(OII_unetpos(:,t)+OI_unetpos(:,t)),'k',y,OII_unet0(:,t)./(OII_unet0(:,t)+OI_unet0(:,t)),'--k',y,OII_unetneg(:,t)./(OII_unetneg(:,t)+OI_unetneg(:,t)),':k',S(kuro),ostreo(kuro)./100,'kx','LineWidth',2,'MarkerSize',12);axis([33.2 34.4 0 1]);ylabel({'Relative abundance';'of OII'},'FontSize',n,'FontName','Helvetica');xlabel('Salinity','FontSize',n,'FontName','Helvetica');
    set(gca,'FontSize',n,'FontName','Helvetica','XTick',33.2:0.2:34.4)
    
    
elseif ftype == 3,
    figure(1);
    hold on
    hline1 = plot(y, out.null(:,1), '--k','LineWidth',2);
    hline2 = plot(y,out.pos(:,1),'k','LineWidth',2);
    hline3 = plot(y,out.neg(:,1),':k','LineWidth',2);
    hline4 = plot(S(kuro),O(kuro,1),'xk','LineWidth',2,'MarkerSize',12);
    
    plot(y, out.null(:,1), '--b',y,out.pos(:,1),'b',y,out.neg(:,1),':b','LineWidth',2); %plot OI output
    plot(y, out.null(:,2), '--r',y,out.pos(:,2),'r',y,out.neg(:,2),':r', 'LineWidth',2); %plot OII output
    plot(S(kuro),O(kuro,1),'xb',S(kuro),O(kuro,2),'xr','LineWidth',2,'MarkerSize',12); % plot observations
    
    ylabel('Abundance (copies ml^{-1})','FontSize',n,'FontName','Helvetica');
    xlabel('Salinity','FontSize',n);set(gca,'FontSize',n)
    set(gca,'FontSize',n,'FontName','Helvetica','XTick',33.2:0.2:34.4, 'LineWidth',2);
    axis([33.2 34.4 0 12000])
    set(gcf, 'Color', 'w');
    
    legend([hline1, hline2, hline3, hline4],'\mu_{NET}=0','\mu_{NET}>0','\mu_{NET}<0','Observations','location','north');%legend 
    legend('boxoff');
    box on
    hold off
    
    elseif ftype == 4,
    kuro_ratio = O(1:31,2)./(O(1:31,1)+O(1:31,2));
    null_ratio = out.null(:,2)/(out.null(:,1)+out.null(:,2));
    pos_ratio = out.pos(:,2)/(out.pos(:,1)+out.pos(:,2));
    neg_ratio = out.neg(:,2)/(out.neg(:,1)+out.neg(:,2));
      
    figure(1);
    hold on
    hline1 = plot(y, null_ratio, '--k','LineWidth',2);
    hline2 = plot(y,pos_ratio,'k','LineWidth',2);
    hline3 = plot(y,neg_ratio,':k','LineWidth',2);
    hline4 = plot(S(1:31),kuro_ratio,'xk','LineWidth',2,'MarkerSize',12);
    
    plot(y, null_ratio, '--k',y,pos_ratio,'k',y,neg_ratio,':k','LineWidth',2); %plot OI output
    plot(S(1:31),kuro_ratio,'xk','LineWidth',2,'MarkerSize',12); % plot observations
    
    ylabel('Relative abundance of OII','FontSize',14,'FontName','Helvetica');
    xlabel('Salinity','FontSize',14);set(gca,'FontSize',14)
    set(gca,'FontSize',14,'FontName','Helvetica','XTick',33.2:0.2:34.4, 'LineWidth',2);
    axis([33.2 34.4 0 1.1])
    set(gcf, 'Color', 'w');
    
    legend([hline1, hline2, hline3, hline4],'\mu_{NET}=0','\mu_{NET}>0','\mu_{NET}<0','Observations','location','north');%legend 
    legend('boxoff');
    box on
    hold off
    
end
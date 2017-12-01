function [ ] = cvEmulPlot( q, eumlPred )

Date=datestr(now,31);


iPlot = logical(prod(eumlPred > 0, 2));
Y_real = eumlPred(iPlot,1:3);
Y_pred = eumlPred(iPlot,4:6);
Y_err = eumlPred(iPlot,7:9);

if (q ~= 3)
    disp('The plots are by default generated for max 3 outputs. Please change the MATLAB code "cvEmulPlot" according to your data.')
    return
end



for i=1:q
    figure
    set(gcf,'PaperUnits','inches');
    width = 3.25;
    height = width*.75;
    set(gcf,'PaperSize',[width height])
    left = 0;
    bottom = 0;
    set(gcf, 'PaperPosition', [left bottom width height]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gca,'LooseInset',[0 0 0.005 0]);
    
    if i==1
        scale = 1e6;
        title_str = 'Melt Pool depth';
        unit_str = '\mum';
        name_str = 'y_depth';
    elseif i==2
        scale = 1e6;
        title_str = 'Melt Pool width';
        unit_str = '\mum';
        name_str = 'y_width';
    else
        scale = 1;
        title_str = 'Peak Temperature';
        unit_str = 'K';
        name_str = 'y_temp';
    end
    
    errorbar(Y_real(:,i)*scale,Y_pred(:,i)*scale,Y_err(:,i)*scale,'s')%,'Color',99*[1,1,1]/255)
    h = refline(1,0);
    set(h,'linewidth',1.5,'color','r')
    grid
    
    if i==1
        ylim([0 60])
    end
    title(['Cross Validation - ',title_str],'FontSize',8.5,'FontName','Times New Roman')
    set(gca,'FontSize',7,'FontName','Times New Roman')
    xlabel(['FEA Simulation [',unit_str,']'],'FontSize',8,'FontName','Times New Roman')
    ylabel(['MVGP Prediction [',unit_str,']'],'FontSize',8,'FontName','Times New Roman')
    
    graphName=[ '.\CodeOutput\' Date(1:10) '_' Date(12:13) Date(15:16)  Date(18:19) '_CV_emulator_'  ];
    
    print([graphName,name_str],'-f1','-dpdf','-r600') % '-dpng' % '-depsc'
    close;
end

end


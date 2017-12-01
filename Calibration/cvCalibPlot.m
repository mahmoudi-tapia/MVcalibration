function [  ] = cvCalibPlot( q, calibPred )

Date=datestr(now,31);


iPlot = logical(prod(calibPred > 0, 2));
yP_RAW = calibPred(iPlot,1:3);
E_yP_actual = calibPred(iPlot,4:6);
V_yP_actual = calibPred(iPlot,7:9);


yP_RAW = calibPred(:,1:3);
E_yP_actual = calibPred(:,4:6);
V_yP_actual = calibPred(:,7:9);



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
    
    if i==2
        scale = 1e6;
        title_str = 'Melt Pool depth';
        unit_str = '\mum';
        name_str = 'y_depth';
    elseif i==1
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
    
    errorbar(yP_RAW(:,i)*scale,E_yP_actual(:,i)*scale,V_yP_actual(:,i)*scale,'s')%,'Color',99*[1,1,1]/255)
    h = refline(1,0);
    set(h,'linewidth',1.5,'color','r')
    grid
    
    % if i==1
    %     ylim([0 60])
    % end
    title(['Validation - ',title_str],'FontSize',8.5,'FontName','Times New Roman')
    set(gca,'FontSize',7,'FontName','Times New Roman')
    xlabel(['Experimental Observation [',unit_str,']'],'FontSize',8,'FontName','Times New Roman')
    ylabel(['MVGP Prediction [',unit_str,']'],'FontSize',8,'FontName','Times New Roman')
    
    graphName=[ '.\CodeOutput\' Date(1:10) '_' Date(12:13) Date(15:16)  Date(18:19) '_CV_calibrator_'  ];
    print([graphName,name_str],'-f1','-dpdf','-r600') % '-dpng' % '-depsc'
    close;
end

end


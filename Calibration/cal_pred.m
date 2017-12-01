%% Prediction
% goes inside calibration directory
clc,clear,close all
% load 2017-09-17_1006_workspace.mat
load 2017-09-17_0243_workspace.mat
% cd 'Calibration code pack'/

% Extracting the prediction data
predData = fullData(validIndex,:);
S = size(predData,1);
N = size(exprData,1);

% Normalize the data  

xP_RAW = predData(:,1:kappa);    %input
% Normalize the input -- between 0 and 1
temp = xS_RAW(:,1:kappa);   % Normalize according to the simulation data
xP = (xP_RAW-repmat(min(temp),length(xP_RAW),1))...
            ./repmat([max(temp)-min(temp)],length(xP_RAW),1);

yP_RAW = predData(:,kappa+1:kappa+q);    %output
% Normalize output -- mean 0
yP = (yP_RAW-repmat(mean(D_RAW),length(yP_RAW),1))...
            ./repmat(std(D_RAW),length(yP_RAW),1); 
        
% Computing terms neeeded for predicting E[yP|yE]
[muDSxP, T_xP] = muDstar2( xP, thetaST, BetaHAT_GLS, D, A, H, xS, r_hat);
[muDSxE, T_xE] = muDstar2( xE, thetaST, BetaHAT_GLS, D, A, H, xS, r_hat);

% Make the vectors stacked
muDSxP = reshape(muDSxP', [S*q,1]);
muDSxE = reshape(muDSxE', [N*q,1]);
yEstacked = reshape(yE', [N*q,1]);

% Compute SigmaEE
C_em = cDstar2( xE, thetaST, T_xE, A, H, r_hat);
SigmaGLS = SigmaHAT2( A, H, D, BetaHAT_GLS);
covI = kron(C_em, SigmaGLS);    % Term I for SigmaEE

C_delta = corrMat2( rDeltaST, xE); 
covII = kron(C_delta, diag(psiDeltaST));    % Term II for SigmaEE

covIII = kron(eye(N), diag(psiEpsST));  % Term III for SigmaEE

SigmaEE = covI + covII + covIII;

% Compute SigmaPE
[ C_em_PE, C_delta_PE  ] = corrPE( xP, xE, thetaST, T_xP, T_xE, A, H, r_hat, rDeltaST);
[ C_em_PP, C_delta_PP  ] = corrPE( xP, xP, thetaST, T_xP, T_xP, A, H, r_hat, rDeltaST);
SigmaPE = kron(C_em_PE,SigmaGLS) + kron(C_delta_PE,diag(psiDeltaST));
SigmaPP = kron(C_em_PP,SigmaGLS) + kron(C_delta_PP,diag(psiDeltaST)) ...
    + kron(eye(S),diag(psiEpsST));

% Kriging for prediction
temp = SigmaEE\(yEstacked - muDSxE);
E_yP = muDSxP + SigmaPE * temp;
E_yP = reshape(E_yP,[q,S])';

temp2 = SigmaEE\SigmaPE';
V_yP = SigmaPP - SigmaPE * temp2;

% -------------------------------
% Plotting the predictions
figure
nplot=q;


E_yP_actual = E_yP .* repmat(std(D_RAW),length(yP_RAW),1) + repmat(mean(D_RAW),length(yP_RAW),1);
V_yP_actual = reshape(diag(V_yP),q,S)' .* repmat(std(D_RAW),length(yP_RAW),1);

% Normal values
% for i=1:q    
%     subplot(nplot,1,i);
%     scatter(1:S, yP(:,i),'filled');
%     grid;  
%     hold on
%     scatter(1:S, E_yP(:,i),'LineWidth',1.5,'MarkerEdgeColor','r');
%     xlabel('Index #')
%     ylabel(['\theta_' num2str(i) ])
% end

% Actual values
for i=1:q    
    subplot(nplot,1,i);
    
    if i ~= 3
        scatter(1:S, 10^6*yP_RAW(:,i),'filled');
        grid; hold on
        scatter(1:S, 10^6*E_yP_actual(:,i),'LineWidth',1.5,'MarkerEdgeColor','r');
    else
        scatter(1:S, yP_RAW(:,i)-273,'filled');
        grid; hold on
        scatter(1:S, E_yP_actual(:,i)-273,'LineWidth',1.5,'MarkerEdgeColor','r');
    end
    
    
    xlabel('Index #')
    
    switch i
        case 1
            ylabel(['\theta_' num2str(i) ': Melt pool width (\mum)' ]);
        case 2
            ylabel(['\theta_' num2str(i) ': Melt pool depth (\mum)' ]);
        case 3
            ylabel(['\theta_' num2str(i) ': Melt pool peak temperature (\circC)' ]);
    end
    
    legend('Experiments','Predictions','Location','Best')
end



% Set the width and height
x0=500;
y0=100;
width=400;
height=1200;
set(gcf,'units','points','position',[x0,y0,width,height])

% Title for the whole plot
set(gcf,'NextPlot','add');
axes;
h = title(['Predictions -- Job ID: ' Date(1:10) '-' Date(12:13) Date(15:16) ]);
set(gca,'Visible','off');
set(h,'Visible','on');

save('cal_pred','yP_RAW','E_yP_actual','V_yP_actual')
% Save the plot file
% graphName=[ '.\CodeOutput\' Date(1:10) '_' Date(12:13) Date(15:16) '_Pred'  ];
% print ('-dtiff', graphName);
% close;

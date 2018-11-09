
function [ thetaSample, thetaMODE, rDeltaST, psiDeltaST, psiEpsST ] = calibratorAfterCV( Date, startTime, r_hat, simData, fullCalibData, P, q, t, MCMC )


iter = MCMC;   % Number of iterations for MCMC
PARAM2 = 'Truncated Normal - Normal proposal for THETA';    % Note on the proposal distributions
PARAM3 = 'Sim. data size: 130, Calib. data 1 with Tp';  % Note on the experimental data
COMMENT1 = 'Experimental data size 24, with 4 validation chose RANDOMLY ';
COMMENT2 = 'Predictions based on POSTERIOR MEAN';



% Load r_hat from the emulator previously built
% load emulatorResult_train60__2017-09-15_1532_workspace__rMode
load emulatorResult_train130__2017-09-15_1411_workspace__rMode

% Acquire the experimental and prediction data
Nfull = size(fullCalibData,1);

% MCMC parameters
thin = 5;
burnin = iter/4;
printPar = 100;     % Show the acceptance rations every X iterations





% --------------------------
%       SIMULATION  DATA
% --------------------------
xS_RAW = simData(:,1:P);    % Collect the input
% Normalize the input -- between 0 and 1
xS = (xS_RAW-repmat(min(xS_RAW),length(xS_RAW),1))...
            ./repmat([max(xS_RAW)-min(xS_RAW)],length(xS_RAW),1);
  
D_RAW = simData(:,P+1:P+q);    % Collect the output matrix
% Normalize output -- mean 0
D = (D_RAW-repmat(mean(D_RAW),length(D_RAW),1))...
            ./repmat(std(D_RAW),length(D_RAW),1);    


% Generate the design matrix H, the correlation matrix A, and the GLS regression coefficients BetaHAT
% We need the results of the emulation (r_hat) for this step
H = designMat (xS, P+1);
A = corrMat2(r_hat, xS);        
BetaHAT_GLS = BetaHAT2( A, H, D);

% -----------------------------
%       EXPERIMENTAL  DATA
% -----------------------------

% Extracting the calibration data
exprData = fullCalibData;
% exprData(validIndex,:) = [];      

% Collect the input
xE_RAW = exprData(:,1:t);    
% Normalize the input -- between 0 and 1
temp = xS_RAW(:,1:t);   % Normalize according to the simulation data
xE = (xE_RAW-repmat(min(temp),length(xE_RAW),1))...
            ./repmat([max(temp)-min(temp)],length(xE_RAW),1);

% Collect the output matrix
yE_RAW = exprData(:,t+1:t+q);    
% Normalize output -- mean 0
yE = (yE_RAW-repmat(mean(D_RAW),length(yE_RAW),1))...
            ./repmat(std(D_RAW),length(yE_RAW),1);    
        
        
%% ------------------------------------------------------   
%       Single Component Metropolis-Hastings (SCMH)
% -------------------------------------------------------

% For testing
% thin = 1;
% burnin = 10;

% Starting values
thetaStart = .5*ones(P-t,1);
rDeltaStart = ones(t,1);
psiDeltaStart = ones(q,1);
psiEpsStart = rand(q,1);

% MCMC sample
thetaSample = zeros((iter-burnin)/thin,P-t);
rDeltaSample = zeros((iter-burnin)/thin,t);
psiDeltaSample = zeros((iter-burnin)/thin,q); 
psiEpsSample = zeros((iter-burnin)/thin,q); 

% Number of acceptances
thetaAcpt = zeros(size(thetaStart));
rDeltaAcpt = zeros(size(rDeltaStart)); 
psiDeltaAcpt = zeros(size(psiDeltaStart));
psiEpsAcpt = zeros(size(psiEpsStart));

% MCMC chain variables
theta = thetaStart;
rDelta = rDeltaStart;
psiDelta = psiDeltaStart;
psiEps = psiEpsStart;


% Set the parameters for the proposal distributions and run the MC chain
propVar = 0.5;

rng default % for reproducability
for i=1:iter
    
    
    
    % -----------------------
    %       theta    
    % -----------------------
    % Update components of theta
    
    for j=1:length(theta)
        
        % 0 -- Adjust the variance for the proposal distibutions
        switch j
            case 1
                propVarTheta=0.2;
            case 2
                propVarTheta=0.3;
            case 3
                propVarTheta=0.1;
        end
            
        
        % 1 -- Compute the current posterior probability        
        postPHI = - post_PHI2(theta, rDelta, psiDelta, psiEps, xE, yE, r_hat, xS, D, A, H, BetaHAT_GLS);

%         ****************************************************************
%                               NORMAL PROPOSAL
%         ****************************************************************
        % 2 -- Update the jth component of theta using the proposal sampler
        newComponent = propNormRnd(theta(j), propVarTheta, 0, 1);     % Normal
        thetaNew = theta;
        thetaNew(j) = newComponent;
        
        % 3 -- Compute proposal probabilities
        prop = -log(propNormPDF(theta, thetaNew, propVarTheta, 0, 1));   % Normal
        propNew = -log(propNormPDF(thetaNew, theta, propVarTheta, 0, 1));      % Normal
%         ****************************************************************

        

%         ****************************************************************
%                               UNIFORM PROPOSAL
%         ****************************************************************
%         2 -- Update the jth component of theta using the proposal sampler
%         newComponent = unifrnd(theta(j)-propVarTheta, theta(j)+propVarTheta);     % UNIFORM
%         thetaNew = theta;
%         thetaNew(j) = newComponent;
%         
%         3 -- Compute proposal probabilities
%         prop = -log(unifpdf(theta, thetaNew-propVarTheta, thetaNew+propVarTheta));   % UNIFORM
%         propNew = -log(unifpdf(thetaNew, theta-propVarTheta, theta+propVarTheta));      % UNIFORM
%         ****************************************************************
        
        % 4 -- Compute new posterior probability
        postPHInew = - post_PHI2(thetaNew, rDelta, psiDelta, psiEps, xE, yE, r_hat, xS, D, A, H, BetaHAT_GLS);
        
        % 5 -- Compute acceptance ratio
        accRatio = exp(-postPHInew - prop + postPHI + propNew);
        
        % 6 -- accept or reject
        if unifrnd(0,1) < accRatio
            theta = thetaNew;
            thetaAcpt(j) = thetaAcpt(j) + 1;
        end
    end
    
    % Store the new sample 
    if mod(i,thin) == 0 && i > burnin
        thetaSample((i-burnin)/thin,:) = theta;
    end

    
    
    % -----------------------
    %       r_delta    
    % -----------------------
    % Update components of r_delta
    for j=1:length(rDelta)
        
        % 1 -- Compute the current posterior probability
        postPHI = - post_PHI2(theta, rDelta, psiDelta, psiEps, xE, yE, r_hat, xS, D, A, H, BetaHAT_GLS);
        
        % 2 -- Update the jth component of r_delta using the proposal sampler
        newComponent = propNormRnd(rDelta(j), propVar, 0, Inf);     % Normal
        rDeltaNew = rDelta;
        rDeltaNew(j) = newComponent;
        
        % 3 -- Compute proposal probabilities
        prop = -log(propNormPDF(rDelta, rDeltaNew, propVar, 0, Inf));   % Normal
        propNew = -log(propNormPDF(rDeltaNew, rDelta, propVar, 0, Inf));      % Normal
        
        % 4 -- Compute new posterior probability
        postPHInew = - post_PHI2(theta, rDeltaNew, psiDelta, psiEps, xE, yE, r_hat, xS, D, A, H, BetaHAT_GLS);
        
        % 5 -- Compute acceptance ratio
        accRatio = exp(-postPHInew - prop + postPHI + propNew);
        
        % 6 -- accept or reject
        if unifrnd(0,1) < accRatio
            rDelta = rDeltaNew;
            rDeltaAcpt(j) = rDeltaAcpt(j) + 1;
        end
        
    end
    
    % Store the new sample 
    if mod(i,thin) == 0 && i > burnin
        rDeltaSample((i-burnin)/thin,:) = rDelta;
    end    
    
    
    
    % -----------------------
    %       psi_delta    
    % -----------------------
    % Update components of psi_delta
    for j=1:length(psiDelta)
        
        % 1 -- Compute the current posterior probability
        postPHI = - post_PHI2(theta, rDelta, psiDelta, psiEps, xE, yE, r_hat, xS, D, A, H, BetaHAT_GLS);
        
        % 2 -- Update the jth component of r_delta using the proposal sampler
        newComponent = propNormRnd(psiDelta(j), propVar, 0, Inf);     % Normal
        psiDeltaNew = psiDelta;
        psiDeltaNew(j) = newComponent;
        
        % 3 -- Compute proposal probabilities
        prop = -log(propNormPDF(psiDelta, psiDeltaNew, propVar, 0, Inf));   % Normal
        propNew = -log(propNormPDF(psiDeltaNew, psiDelta, propVar, 0, Inf));      % Normal
        
        % 4 -- Compute new posterior probability
        postPHInew = - post_PHI2(theta, rDelta, psiDeltaNew, psiEps, xE, yE, r_hat, xS, D, A, H, BetaHAT_GLS);
        
        % 5 -- Compute acceptance ratio
        accRatio = exp(-postPHInew - prop + postPHI + propNew);
        
        % 6 -- accept or reject
        if unifrnd(0,1) < accRatio
            psiDelta = psiDeltaNew;
            psiDeltaAcpt(j) = psiDeltaAcpt(j) + 1;
        end
        
    end
    
    % Store the new sample 
    if mod(i,thin) == 0 && i > burnin
        psiDeltaSample((i-burnin)/thin,:) = psiDelta;
    end 
    
    
    

    
    
    % -----------------------
    %       psi_eps    
    % -----------------------
    % Update components of psi_eps
    for j=1:length(psiEps)
        
        % 0 -- Adjust the variance for the proposal distibutions
        if j==2
            propVarEps = .05;
        else
            propVarEps = 0.5;
        end
            
        
        
        % 1 -- Compute the current posterior probability
        postPHI = - post_PHI2(theta, rDelta, psiDelta, psiEps, xE, yE, r_hat, xS, D, A, H, BetaHAT_GLS);
        
        % 2 -- Update the jth component of r_delta using the proposal sampler
        newComponent = propNormRnd(psiEps(j), propVarEps, 0, Inf);     % Normal
        psiEpsNew = psiEps;
        psiEpsNew(j) = newComponent;
        
        % 3 -- Compute proposal probabilities
        prop = -log(propNormPDF(psiEps, psiEpsNew, propVarEps, 0, Inf));   % Normal
        propNew = -log(propNormPDF(psiEpsNew, psiEps, propVarEps, 0, Inf));      % Normal
        
        % 4 -- Compute new posterior probability
        postPHInew = - post_PHI2(theta, rDelta, psiDelta, psiEpsNew, xE, yE, r_hat, xS, D, A, H, BetaHAT_GLS);
        
        % 5 -- Compute acceptance ratio
        accRatio = exp(-postPHInew - prop + postPHI + propNew);
        
        % 6 -- accept or reject
        if unifrnd(0,1) < accRatio
            psiEps = psiEpsNew;
            psiEpsAcpt(j) = psiEpsAcpt(j) + 1;
        end
        
    end
    
    % Store the new sample 
    if mod(i,thin) == 0 && i > burnin
        psiEpsSample((i-burnin)/thin,:) = psiEps;
    end   
    

    
    % Display the iteration number and acceptance ratios
    if mod(i,printPar) == 0
        % display current samples
        fprintf('Iteration: %d\n', i)
        fprintf('Acceptance ratios for theta: %.2f %.2f %.2f \n', thetaAcpt/i)
        fprintf('Acceptance ratios for r_delta: %.2f %.2f \n', rDeltaAcpt/i)
        
        if q==3
        fprintf('Acceptance ratios for psi_delta: %.2f %.2f %.2f \n', psiDeltaAcpt/i)
        fprintf('Acceptance ratios for psi_eps: %.2f %.2f %.2f \n', psiEpsAcpt/i)
        end
       
        if q==2
        fprintf('Acceptance ratios for psi_delta: %.2f %.2f \n', psiDeltaAcpt/i)
        fprintf('Acceptance ratios for psi_eps: %.2f %.2f \n', psiEpsAcpt/i)
        end
        
    end
    
end

% Collecting the calibrated parameters
% thetaST = mean(thetaSample)';
rDeltaST = mean(rDeltaSample)';
psiDeltaST = mean(psiDeltaSample)';
psiEpsST = mean(psiEpsSample)';


%% TRACE PLOTS & HISTOGRAMS

% ------------------------------------
%       Trace plots for all parameters
% ------------------------------------
figure
nplot=max([P-t,t,q]);

for i=1:P-t    
    subplot(4,nplot,i);
    plot(thetaSample(:,i));
    grid;  
    xlabel('Index #');
    ylabel(['\theta_' num2str(i)])
end

for i=nplot+1:nplot+t
    subplot(4,nplot,i);
    plot(rDeltaSample(:,i-nplot));
    grid;  
    xlabel('Index #');
    ylabel(['r_{\delta' num2str(i-nplot) '}'])
end

for i=2*nplot+1:2*nplot+q    
    subplot(4,nplot,i);
    plot(psiDeltaSample(:,i-2*nplot));
    grid;  
    xlabel('Index #');
    ylabel(['\psi_{\delta' num2str(i-2*nplot) '}'])
end


for i=3*nplot+1:3*nplot+q    
    subplot(4,nplot,i);
    plot(psiEpsSample(:,i-3*nplot));
    grid;  
    xlabel('Index #');
    ylabel(['\psi_{\epsilon' num2str(i-3*nplot) '}'])
end

x0=50;
y0=100;
width=1200;
height=900;
set(gcf,'units','points','position',[x0,y0,width,height])

% Title for the whole plot
set(gcf,'NextPlot','add');
axes;
h = title('Trace plots');
set(gca,'Visible','off');
set(h,'Visible','on');

% Write info on the plot
tempX = xlim;
tempY = ylim;
xText = tempX(2)*9/12;
yText = tempY(2)*11/16;
strLine1 = ['Job ID: ' Date(1:10) '-' Date(12:13) Date(15:16) ];
strLine2 = ['Number of MCMC iterations: ' num2str(iter)];
strLine3 = ['Experimental data: ' PARAM3 ];
strLine4 = ['Proposal: ' PARAM2 ];
strLine5 = ['Comment 1: ' COMMENT1];
strLine6 = ['Comment 2: ' COMMENT2];

text(xText,yText, strLine1);
text(xText,yText*.95, strLine2);
text(xText,yText*.90, strLine3);
text(xText,yText*.85, strLine4);
text(xText,yText*.80, strLine5);
text(xText,yText*.75, strLine6);


% Save the plot file
graphName=[ '.\CodeOutput\' Date(1:10) '_' Date(12:13) Date(15:16)  '_TracePlots'  ];
print ('-dtiff', graphName);
close;


% ------------------------------------
%       Histograms
% ------------------------------------
figure
nplot=max([P-t,t,q]);

nHist=100;
for i=1:P-t    
    subplot(4,nplot,i);
    hist(thetaSample(:,i),nHist);
    grid;  
    xlabel(['\theta_' num2str(i) ' (mean: ' num2str(mean(thetaSample(:,i))) ')'])
end

for i=nplot+1:nplot+t
    subplot(4,nplot,i);
    hist(rDeltaSample(:,i-nplot),nHist);
    grid;  
    xlabel(['r_{\delta' num2str(i-nplot) '}' ' (mean: ' num2str(mean(rDeltaSample(:,i-nplot))) ')'])
end

for i=2*nplot+1:2*nplot+q    
    subplot(4,nplot,i);
    hist(psiDeltaSample(:,i-2*nplot),nHist);
    grid;
    xlabel(['\psi_{\delta' num2str(i-2*nplot) '}' ' (mean: ' num2str(mean(psiDeltaSample(:,i-2*nplot))) ')'])
end

for i=3*nplot+1:3*nplot+q    
    subplot(4,nplot,i);
    hist(psiEpsSample(:,i-3*nplot),nHist);
    grid;
    xlabel(['\psi_{\epsilon' num2str(i-3*nplot) '}' ' (mean: ' num2str(mean(psiEpsSample(:,i-3*nplot))) ')'])
end

x0=50;
y0=100;
width=1200;
height=900;
set(gcf,'units','points','position',[x0,y0,width,height])

% Title for the whole plot
set(gcf,'NextPlot','add');
axes;
h = title('Histograms');
set(gca,'Visible','off');
set(h,'Visible','on');

% Write info on the plot
tempX = xlim;
tempY = ylim;
xText = tempX(2)*9/12;
yText = tempY(2)*11/16;
strLine1 = ['Job ID: ' Date(1:10) '-' Date(12:13) Date(15:16) ];
strLine2 = ['Number of MCMC iterations: ' num2str(iter)];
strLine3 = ['Experimental data: ' PARAM3 ];
strLine4 = ['Proposal: ' PARAM2 ];
strLine5 = ['Comment 1: ' COMMENT1];
strLine6 = ['Comment 2: ' COMMENT2];

text(xText,yText, strLine1);
text(xText,yText*.95, strLine2);
text(xText,yText*.90, strLine3);
text(xText,yText*.85, strLine4);
text(xText,yText*.80, strLine5);
text(xText,yText*.75, strLine6);

% Save the plot file
graphName=[ '.\CodeOutput\' Date(1:10) '_' Date(12:13) Date(15:16) '_Hist'  ];
print ('-dtiff', graphName);
close;

%% Plot only theta
THETA = thetaSample;

minVal = [.2, .4, 1];
maxVal =[.7, .9, 25];

THETA = bsxfun(@times,THETA,maxVal-minVal)+minVal;

% Mean, mode, and std. dev
% thetaMEAN = mean(THETA)';
thetaMODE = zeros(P-t,1);
% thetaSD = std(THETA)';

% bin = [20,20,40];
bin = [60,60,60];
npoints = size(THETA,1);
l = minVal-minVal*.075;
u = maxVal+maxVal*.05;


% Subplot
figure
% set(gcf,'PaperUnits','inches');
% width = 3.25;
% height = 3*width;
% set(gcf,'PaperSize',[width height])

x0=50;
y0=100;
width=400;
height=1200;
set(gcf,'units','points','position',[x0,y0,width,height])

%----------------
% SUBPLOT
%----------------
for i=1:q
  

    
        subplot(3,1,i)
    
    x = THETA(:,i);
    [f,xi] = ksdensity(x,'npoints',npoints,'support',[l(i),u(i)]);
    y = histcounts(x,bin(i),'normalization','pdf');
    h = histogram(x,bin(i),'EdgeColor','w','normalization','pdf');
    %     h.FaceColor = 82*[1,1,1]/255; %black and white
    hold on
    plot(xi,f,'Color','r','LineWidth',2)
    grid
    
    
    [~,modeIndex]=max(f);
    thetaMODE(i)=xi(modeIndex);
    
    xlim([minVal(i),maxVal(i)])
    if i==1
        xlim([minVal(i),maxVal(i)])
        strunit = '%';
        strname = 'Porosity';
    elseif i==2
        xlim([.4,.9])
        strunit = '%';
        strname = 'Absorptivity';
    else
%         xlim([3,11])
        xlim([1,25])
        strunit = 'dimensionless';
        strname = 'Liquid Thermal Conductivity Coef.';
    end
    xlabel(['\theta_{',num2str(i),'} [',strunit,']'],'FontSize',8,'FontName','Times New Roman')
    title(strname,'FontSize',8.5,'FontName','Times New Roman')
    ylabel('Probability','FontSize',8,'FontName','Times New Roman')
    set(gca,'FontSize',7,'FontName','Times New Roman')
    
end
%     set(gcf,'color','none')
%     print('thetaPost','-f1','-dpdf','-r600')
%     print ('-dpdf', 'thetaPostAll3', '-r600');
graphName=[ '.\CodeOutput\' Date(1:10) '_' Date(12:13) Date(15:16) '_thetaPostAll'  ];
print ('-dpdf', graphName, '-r600');
close;




%% Save the workspace

% ---------------------    
%       QUICK INFO
% ---------------------
f1 = 'NoInput';
f2 = 'NoOutput';
f3 = 'Data';
f4 = 'max_iteration';
f5 = 'Proposal';

% Save important info of the code
f1val = P;
f2val = q;
f3val = PARAM3;
f4val = iter;
f5val = PARAM2;

simSummary = struct(f1, f1val,...
                    f2, f2val,...
                    f3, f3val,...
                    f4, f4val,...
                    f5, f5val);

elapsedTime = toc(startTime);

Name1=[ '.\CodeOutput\' Date(1:10) '_' Date(12:13) Date(15:16)  '_workspace' ];
save(Name1)
end



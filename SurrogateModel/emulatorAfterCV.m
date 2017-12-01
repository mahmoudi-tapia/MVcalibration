%% Last update: 09/13/2017

% PROJECT NAME: Multi-variate UQ
% SUBPROJECT:   Multi-output emulator

% This script generates the emulator
% Please refer to the document "Implementation" to understand the procedure
    
function [ rMode, rSample ] = emulatorAfterCV( trainingData, p, q, PARAM1, PARAM2, PARAM3 )


Date=datestr(now,31);
tic




% Proceed with the reduced data
inputRAW = trainingData(:,1:p);    % Collect the input
% Normalize the input -- between 0 and 1
inputNORM = bsxfun(@rdivide,bsxfun(@minus,inputRAW,min(inputRAW)),max(inputRAW)-min(inputRAW));
  
outputRAW = trainingData(:,p+1:p+q);    % Collect the output matrix
% Normalize output -- mean 0
D = bsxfun(@rdivide,bsxfun(@minus,outputRAW,mean(outputRAW)),std(outputRAW));    

% Generate the design matrix
% Here we take m = 1
% H = designMat (inputNORM, 1);
H = designMat (inputNORM, p+1);

%% Verify the posterior dist. function at a random point
% % rTest = [11, 0.5, 2, 1.5, 2]'
% rTest = ones(p,1)
% 
% % Test you get a reasonable correlation matrix without too many zero
% % elements
% [A_test,A_test_chol] = corrMat2(rTest, inputNORM);
% temp1=reshape(A_test, length(trainingData)^2,1);
% temp2= temp1>0.05;
% fprintf('Number of nonzero elements in A = %i \n', sum(temp2))
% 
% % Test the posterior probability
% pi_test = log_pi_r_D2( rTest, D, A_test, A_test_chol, H, q )
% 
% % Test the mu double star
% [ BetaHAT_test ] = BetaHAT2( A_test_chol, H, D);
% xTest = rand(10,p);
% [ mu_test ] = muDstar2( xTest, BetaHAT_test, D, A_test_chol, H, inputNORM, rTest)


%% Single Component Metropolis-Hastings (SCMH)
% Use MH algorithm to generate a sample from the posterior distribution
% of vector r which is the vector of positive roughness parameters

rng shuffle %default % for reproducability % ***** Gustavo: removed default


% Starting vector
rStart = ones(p,1);
% rStart = rTest;

% MCMC iterations and parameters
iter = PARAM1;
thin = 5;
burnin = iter/4;

% MCMC sample
rSample = zeros((iter-burnin)/thin,p);
rAccept = zeros(1, p); % number of acceptances

% MCMC chain
r = rStart;
propVar = 0.5;
% Compute cov matrix posterior probability
[A,Achol] = corrMat2(r, inputNORM);
rPost = -(log_pi_r_D2(r, D, A, Achol, H, q ));
for i=1:iter

    % Update every component of the vector r
    for j=1:length(r)
        
        % Update the jth component from the proposal sampler
%         newComponent = propExpRnd(r(j));   % Exponential
        newComponent = propNormRnd(r(j), propVar, 0, Inf);     % Normal
        rNew = r;
        rNew(j) = newComponent;
        
        % Compute proposal probabilities
%         prop = -log(propExpPDF(r, rNew));   % Exponential
        prop = -log(propNormPDF(r, rNew, propVar, 0, Inf));   % Normal
%         propNew = -log(propExpPDF(rNew, r));                    % Exponential
        propNew = -log(propNormPDF(rNew, r, propVar, 0, Inf));      % Normal
        
        % Compute new cov matrix posterior probability
        [Anew, Acholnew] = corrMat2(rNew, inputNORM); % **** Gustavo: a typo, using inputRAW insted of inputNORM
        %rPostNew = -log(pi_r_D(rNew, D, A, H, q ));
        rPostNew = -(log_pi_r_D2(rNew, D, Anew, Acholnew, H, q )); %*** Gustavo: a typo, using A insted of Anew
        
        % Compute acceptance ratio
        accRatio = exp(-rPostNew - prop + rPost + propNew);
        
        % accept or reject
        if unifrnd(0,1) < accRatio
            r = rNew;
            rAccept(j) = rAccept(j) + 1;
            rPost = rPostNew;
        end
        
    end
    
    % Store the new sample 
    if mod(i,thin) == 0 && i > burnin
        rSample((i-burnin)/thin,:) = r;
    end
    
    if mod(i,100) == 0
        % display current samples
        fprintf('Iteration: %d\n', i)
        fprintf('Acceptance ratios: %.2f %.2f %.2f %.2f %.2f \n', rAccept/i)
    end
    
end


elapsedTime = toc;

%% Plotting the histograms and kernel densities

% Declare the modes for posteriors
rMode = zeros(p,1);


bin = 30;
npoints = size(rSample,1);
% l = minVal-minVal*.05;
% u = maxVal+maxVal*.02;

figure
for i=1:p
   
    
%     set(gcf,'PaperUnits','inches');
%     width = 3;
%     height = width*.75;
%     set(gcf,'PaperSize',[width height])
%     left = 0;
%     bottom = 0;
%     set(gcf, 'PaperPosition', [left bottom width height]);
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gca,'LooseInset',[0 0 0.005 0]);
    
    subplot(p,1,i)
    
    x = rSample(:,i);
    % [f,xi] = ksdensity(x,'npoints',npoints,'support',[l(i),u(i)]);
    [f,xi] = ksdensity(x,'npoints',npoints);
    y = histcounts(x,bin,'normalization','pdf');
    h = histogram(x,bin,'EdgeColor','w','normalization','pdf');
    h.FaceColor = 82*[1,1,1]/255; %black and white
    hold on
    plot(xi,f,'Color','b','LineWidth',2)
    grid
    % xlim([minVal(i),maxVal(i)])
    
    % Find the mode
    [~,modeIndex]=max(f);
    rMode(i)=xi(modeIndex);
    
    xlabel(['r_{',num2str(i),'}' ' -- (mean: ' num2str(mean(rSample(:,i))) ...
        ', mode: ' num2str(rMode(i)) ,')'],...
    'FontSize',8,'FontName','Times New Roman')
%     title(['r_{',num2str(i),'}'],'FontSize',8.5,'FontName','Times New Roman')
    ylabel('Probability','FontSize',8,'FontName','Times New Roman')
    set(gca,'FontSize',7,'FontName','Times New Roman')

    
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
h = title(['Posteriors for r -- Job ID: ' Date(1:10) '-' Date(12:13) Date(15:16) ]);
set(gca,'Visible','off');
set(h,'Visible','on');


% Save the plot file
graphName=[ '.\CodeOutput\' Date(1:10) '_' Date(12:13) Date(15:16)  Date(18:19) '_rHist'  ];
print ('-dtiff', graphName);
close;

%% Save the workspace

% Save important info of the code
f1 = 'NoInput';
f2 = 'NoOutput';
f3 = 'Data';
f4 = 'max_iteration';
f5 = 'Proposal';


% Save important info of the code
f1val = p;
f2val = q;
f3val = PARAM3;
f4val = iter;
f5val = PARAM2;

simSummary = struct(f1, f1val,...
                    f2, f2val,...
                    f3, f3val,...
                    f4, f4val,...
                    f5, f5val);
                
Name1=[ '.\CodeOutput\' Date(1:10) '_' Date(12:13) Date(15:16) Date(18:19)  '_final_emulator_workspace' ];
save(Name1)











end
function draw_fitBOLD(results,fitstats,fit,lambda_BIC)

if ~isempty(strfind(fit.dependence,'fused'))
    if AR~=0
        figure('pos',[10 10 700 1000]);subplot(4,1,1);
        plot(results.y,'b'); hold on; plot(fitstats.lagged_yhat,'r');
        title('Lagged predictions against actual ratings');
        dataNames_1 = {'Actual data','lagged predicted data'};
        legend(dataNames_1);
        ylabel([sprintf('R2 = %.2f',fitstats.lagged_BOLDr2_fused)...
            sprintf('   LSE = %.2f',fitstats.lagged_BOLD_lse)]);
        xlabel([sprintf('t-statistic = %.2f ',fitstats.lagged_BOLD_tstat)...
            sprintf('   p-value = %.2f',fitstats.lagged_BOLD_pval)]);
        subplot(4,1,2);
        plot(results.y,'b'); hold on; plot(fitstats.lagterms_pred,'r');
        title(['Autoregression(' num2str(fitstats.nlag) ') Lag terms predictions against actual ratings']);
        dataNames_1 = {'Actual data','predictions of lag terms'};
        legend(dataNames_1);
        ylabel([sprintf('R2 = %.2f',fitstats.lagterms_r2)...
            sprintf('   LSE = %.2f',fitstats.lagterms_lse)]);
        xlabel([sprintf('t-statistic = %.2f ',fitstats.lagterms_tstat)...
            sprintf('   p-value = %.2f',fitstats.lagterms_pval)]);
        
        subplot(4,1,[3 4])
    end
    figure;
    for g=1:size(fitstats.lambda_BIC,1)
        lambda_BIC(g,:)=[fitstats.lambda_BIC{g,1};];
    end
    cVals = unique(lambda_BIC(:,4));
    for i = 1:numel(cVals)                     % For every one of those unique values
        indices = find(lambda_BIC(:,4) == cVals(i));         % Find the corresponding indices
        scatter3(lambda_BIC(indices,1),lambda_BIC(indices,2),lambda_BIC(indices,3),100,'filled') % Plot
        hold on
    end
    lgd1=cell(1,0);
    for lgd=1:size(cVals,1)
        lgd1=horzcat(lgd1,{['BIC=  ' num2str(round(cVals(lgd,1)),2)]});
    end
    legend(lgd1)
    % plot([fitstats.lambda_BIC(:,1)], [fitstats.lambda_BIC(:,4)])
    xlabel('Fused lasso lambda');  ylabel('L1 norm lambda "Lasso"');zlabel('L2 norm lambda "Ridge"'); title('Flat Fused Lasso(FISTA optimization): BIC for each Regularization parameter Lambda');
    
end
figure('pos',[10 10 700 900]);
subplot(2,1,1)
plot(results.y,'b');
hold on; plot(results.train_yhat,'r')
dataNames = {'Actual data','predicted data'};
if ~isempty(strfind(fit.dependence,'fused'))
    dataNames{1,2}=['AR(' num2str(fitstats.nlag) ') corrected predictions'];
end
legend(dataNames)


if fit.AR==1
    xlabel([sprintf('t-statistic = %.2f ',fitstats.model_tstat)...
        sprintf('   p-value = %.2f',fitstats.model_pval)...
        sprintf('   F-statistic = %.2f ',fitstats.ftest)...
        sprintf('   F test p-val = %.2f',fitstats.fprob)]);
    if fitstats.nlag >0
        txt=['  Sims Corrected p-value: ' num2str(round(fitstats.Sims_dof_correced_probability,3)) ' Granger marginal probability: ' num2str(round(fitstats.fprob,3))];
    else
        txt=[' Granger marginal probability ' num2str(round(fitstats.fprob,3))];
    end
    
    title(['Subject: ' num2str(results.sub) txt], 'interpreter','none')
else
    if ~isempty(strfind(fit.dependence,'step')) || ~isempty(strfind(fit.dependence,'OLS'))
        xlabel([sprintf('F-statistic = %.2f ',results.model_tstat)...
            sprintf('p-value = %.2f',results.model_pval)]);
    else
        xlabel([sprintf('t-statistic = %.2f ',fitstats.model_tstat)...
            sprintf('   p-value = %.2f',fitstats.model_pval)]);
    end
end
if ~isempty(strfind(fit.dependence,'step'))
    ylabel([sprintf('R2 = %.2f',results.rsqr)...
        sprintf('   SSE = %.2f',fitstats.SSE)]);
end

subplot(2,1,2)
if ~isempty(strfind(fit.dependence,'step'))
    Ccl=zeros(1,0);
    grpcl=zeros(1,0);
    CI=coefCI(fitstats);
    try
    for get_box_cl= 1:size(results.ROIs,2)
        Ccl = [Ccl [CI(get_box_cl,:)] ];
        grpcl = [grpcl ((get_box_cl-1)+zeros(1,size(CI,2)))];
    end
    scatter(1:size([results.betas]),results.betas)
    hold on;
    boxplot(Ccl,grpcl);
    title(['Parameters estimates & 95% Confidence interval'], 'interpreter','none');
    end
else
    bar(results.betas)
    title(['Parameters estimates'], 'interpreter','none');
end

lbl=cell(1,size(results.ROIs,2));
for x_lbl=1:size(results.ROIs,2)
    if isempty(strfind(fit.dependence,'step'))
        lbl{1,x_lbl}=[results.ROIs{1,x_lbl} ' at ' results.timings{1,x_lbl}];
    else
        lbl=results.ROIs;
    end
    try
        if results.pval(x_lbl,1) <0.05
            lbl{1,x_lbl}=[lbl{1,x_lbl} '*'];
        end
    end
end
set(gca,'XTick',1:size(results.ROIs,2));
try xlabel_oblique(lbl); end
drawnow

figure('pos',[10 10 700 900]);
subplot(2,1,1)
title('Residuals Autocorrelation', 'interpreter','none');
autocorr_c(fitstats.Residuals.Raw,5);
%dataNames = {'Raw residuals','Standardized residuals'};
%legend(dataNames)
subplot(2,1,2)
histfit(fitstats.Residuals.Raw);
title('Residuals histogram', 'interpreter','none');


fprintf('Saving figures \n')
h=findobj('type','figure'); % find the handles of the opened figures
folder=fit.dest_dir;  % Desination folder

for k=1:numel(h)
    filename=sprintf('%d.jpg',k);
    filename2=sprintf('%d.fig',k);
    fileset1= name_fitBOLD(fit);
    figure_name=[fileset1 '_' filename];
    figure_name2=[fileset1 '_' filename2];
    file=fullfile(folder,figure_name);
    file2=fullfile(folder,figure_name2);
    saveas(h(k),file)
    saveas(h(k),file2)
    close
end
end

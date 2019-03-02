function [rescell,rescell_Bmodel]=draw_fitBOLDp(results,fitstats,fit)
modelName=name_fitBOLD(fit);

if ~isempty(strfind(fit.dependence,'fused'))
    if fit.AR~=0
        model_b=find([results(:).prediction_r2_lagcorct]==max([results(:).prediction_r2_lagcorct]));
    else
        model_b=find([results(:).prediction_r2]==max([results(:).prediction_r2]));
    end
    
    for m = 1:size(results,2)
        figure('pos',[10 10 1200 900]);
        subplot(3,1,1)
        try  plot(results(m).train_y,'b');
            hold on; plot(results(m).train_yhat,'r'); end
        dataNames = {'Actual data','training data predictions'};
        legend(dataNames)
        if fit.AR~=0
            rf=['Model No(' num2str(m) ') Training data OLSAR(' num2str(results(m).nlag) ') regression'];
        else
            rf=['Model No(' num2str(m) ') Training data OLSAR regression'];
        end
        title(rf, 'interpreter','none')
        ylabel(sprintf('Regression R2 = %.2f',results(m).rsqr));
        xlabel(['BIC=' num2str(results(m).bic) '  AIC=' num2str(results(m).aic)])
        subplot(3,1,2)
        plot(results(m).test_y,'b');
        hold on; plot(results(m).test_yhat,'r')
        dataNames = {'Actual data','test data predictions'};
        legend(dataNames)
        ylabel([sprintf('Test prediction R2 = %.2f',results(m).prediction_r2) sprintf('   p-value = %.2f',results(m).prediction_pval)]);
        rf=['Model No(' num2str(m) ') Test data predictions versus actual test data'];
        title(rf, 'interpreter','none')
        subplot(3,1,3)
        bar(results(m).betas(1:size(results(m).ROIs,2)))
        title([ 'Model No(' num2str(m) ') '  modelName ' - parameters'], 'interpreter','none')
        lbl=cell(1,size(results(m).ROIs,2));
        if isempty(strfind(fit.dependence,'step'))
            for x_lbl=1:size(results(m).ROIs,2)
                lbl{1,x_lbl}=[results(m).ROIs{1,x_lbl} ' at ' results(m).timings{1,x_lbl}];
                
                if results(m).pval(x_lbl,1) <0.05
                    lbl{1,x_lbl}=[lbl{1,x_lbl} '*'];
                end
            end
        else
            for x_lbl=1:size(results(m).ROIs,2)
                if results(m).pval(x_lbl,1) <0.05
                    lbl=results(m).ROIs;
                    lbl{1,x_lbl}=[lbl{1,x_lbl} '*'];
                end
                
            end
        end
        set(gca,'XTick',1:size(results(m).ROIs,2));
        try xlabel_oblique(lbl); end
        
        %*****Plotting Best model cross validations
        figure('Name','Best model - Cross Validations','pos',[10 10 1200 900]);
        subplot(4,1,1)
        plot(fitstats(m).y,'b');
        hold on; plot(fitstats(m).yhat,'r')
        dataNames = {'Actual data','training data predictions'};
        legend(dataNames)
        if fit.AR~=0
            rf=['Best model - Cross Validation No(' num2str(m)  ') Training data OLSAR(' num2str(fitstats(m).nlag) ') regression'];
        else
            rf=['Best model - Cross Validation No(' num2str(m)  ') Training data OLS regression'];
        end
        title(rf, 'interpreter','none')
        ylabel(sprintf('Regression R2 = %.2f',fitstats(m).rsqr));
        xlabel(['BIC=' num2str(fitstats(m).bic) '  AIC=' num2str(fitstats(m).aic)])
        
        subplot(4,1,2)
        plot(fitstats(m).testB_y,'b');
        hold on; plot(fitstats(m).testB_yhat,'r')
        dataNames = {'Actual data','test data predictions'};
        legend(dataNames)
        ylabel([sprintf('prediction R2 = %.2f',fitstats(m).prediction_r2) sprintf('   p-value = %.2f',fitstats(m).prediction_pval)]);
        rf=['Best model - Cross Validation No(' num2str(m) ') Test data predictions versus actual test data'];
        title(rf, 'interpreter','none')
        if fit.AR~=0
            subplot(4,1,3)
            plot((fitstats(m).testB_y-fitstats(m).lag_predB),'b');
            hold on; plot(fitstats(m).testB_yhat,'r')
            dataNames = {'Actual data - lag corrected','test data predictions'};
            legend(dataNames)
            ylabel([sprintf('prediction R2 = %.2f',fitstats(m).prediction_r2_lagcorct) sprintf('   p-value = %.2f',fitstats(m).prediction_pval_lagcorct)]);
            rf=['Best model - Cross Validation No(' num2str(m) '): Test data predictions vs. lag(' num2str(fitstats(m).nlag) ') - corrected actual test data'];
            title(rf, 'interpreter','none')
            subplot(4,1,4)
        else
            subplot(4,1,3)
        end
        bar(fitstats(m).beta(1:size(results(model_b).ROIs,2)))
        title(['Best model - Cross Validation No(' num2str(m) modelName ' - parameters,'], 'interpreter','none')
        lbl=cell(1,size(results(model_b).ROIs,2));
        for x_lbl=1:size(results(model_b).ROIs,2)
            lbl{1,x_lbl}=[results(model_b).ROIs{1,x_lbl} ' at ' results(model_b).timings{1,x_lbl}];
            
            if fitstats(m).pval(x_lbl,1) <0.05
                lbl{1,x_lbl}=[lbl{1,x_lbl} '*'];
            end;
        end
        set(gca,'XTick',1:size(results(model_b).ROIs,2));
        try xlabel_oblique(lbl); end
    end
    [rescell,~,rescell_Bmodel]=orgres_fitBOLD(results,fitstats,fit,0);
elseif ~isempty(strfind(fit.dependence,'step'))
  try
      if fit.regress.validate
   figure('pos',[10 10 1200 1000]); 
   for i=1:size(results.predictedresponse,2)
       subplot(size(results.predictedresponse,2),1,i)
       plot(results.testedresponse{1,i},'b');
       hold on;
       plot(results.predictedresponse{1,i},'r')
       ylabel(['R2=' num2str(results.validationR2(i)) '  RMSE=' num2str(results.validationRMSE(i))])
   end
      end
  end
    figure('pos',[10 10 1200 1000]);
     subplotnum=5;
     nobs=size(results.y,1);
    divvec=0: round(size(results.y,1)/ subplotnum):size(results.y,1);
    divvec(end+1)=nobs;
    for div=1: subplotnum
    yhat=results.BOLDpred(divvec(div)+1:divvec(div+1),1);
    yhatci=results.yhatci(:,divvec(div)+1:divvec(div+1))';
    y=results.y(divvec(div)+1:divvec(div+1),1);
    subplot(subplotnum,1,div)
    plot(y,'b');
    hold on; plot(yhat,'r');
    hold on; plot(yhatci,'k')
    dataNames = {'Actual data','mean predicted data','95% Confidence interval'};
    legend(dataNames)
    rf='Pooled data from all subjects - 95% CI';
    title(rf, 'interpreter','none') 
    ylabel(sprintf('R2 = %.2f',results.rsqr));
    xlabel(['BIC=' num2str(results.bic) '  AIC=' num2str(results.aic)])
    end
    if ~isempty(results.coefci)
    figure('pos',[10 10 1200 1200]);
    %subplotnum2=2;
    nocoef=size(results.coefci,2);
    %subloc=1:2:subplotnum2;
    %divvec2=0: round(nocoef/subplotnum2):nocoef;
    %divvec2(end+1)=nocoef;
    %for div2=1: subplotnum2
    %if div2==(subplotnum2/2)+1
    %figure('pos',[10 10 1200 1200]);
    %end
    Ccl=zeros(1,0);
    grpcl=zeros(1,0);
    %CI=results.coefci{[2,3],(divvec2(div2)+1:divvec2(div2+1))}';
    %sub_lbl={results.predictorROIs{1,(divvec2(div2)+1:divvec2(div2+1))}};
     CI=results.coefci{[2,3],:}';
    sub_lbl={results.predictorROIs{1,:}};
    try
    for get_box_cl= 1:size(CI,1)
        Ccl = [Ccl [CI(get_box_cl,:)] ];
        grpcl = [grpcl ((get_box_cl-1)+zeros(1,size(CI,2)))];
    end
    end
    %if div2>(subplotnum2/2); subloc_i=div2-(subplotnum2/2); else subloc_i=div2;end
    %subplot(subplotnum2,1,subloc(subloc_i))
    boxplot(Ccl,grpcl);
    title(['Parameters estimates & 95% Confidence interval'], 'interpreter','none');
    lbl=sub_lbl;
    set(gca,'XTick',1:size(lbl,2));
    try xlabel_oblique(lbl); end
    %end
    if fit.use_subROIs==1 && results.model_in==1
    coordsunpacked=[];
    for unpack_coords=1:size(results.sub_coord,2)
    coordsunpacked=[coordsunpacked,results.sub_coord{1,unpack_coords}];
    end
    coordsall=horzcat(coordsunpacked,coordsunpacked);
    coordsall_in={coordsall{1,(sum(results.betas,1)~=0)}};
    linbetas=results.coefci{1,results.linearityidx(1,:)==1};
    %linrois=results.predictorROIs{1,results.linearityidx(1,:)==1};
    lincoords={coordsall_in{1,results.linearityidx(1,:)==1}};
    quadbetas=results.coefci{1,results.linearityidx(2,:)==1};
    %quadrois=results.predictorROIs{1,results.linearityidx(2,:)==1};
    quadcoords={coordsall_in{1,results.linearityidx(2,:)==1}};
    imagescoords={lincoords,quadcoords};
    imagesint={linbetas,quadbetas};
    imagelbl={'linear','quadratic'};
        
        cd(fit.dest_dir)
        load('mm3.mat');
        img_mat=0;
        for image_i=1:size(imagescoords,2)
        for subROIs_i= 1:size(imagescoords{1,image_i},2)
            ROIs_array{1,subROIs_i}= maroi_pointlist(struct('XYZ',imagescoords{1,image_i}{1,subROIs_i},'mat',mm3,'roithresh',0,'vals',zeros(1,size(imagescoords{1,image_i}{1,subROIs_i},2))+imagesint{1,image_i}(subROIs_i),'binarize',0),'vox');
            img_mat=img_mat+ROIs_array{1,subROIs_i};
        end
        clear ROIs_array
        modelName=name_fitBOLD(fit);
        lbl_roi=[modelName 'betas_'  imagelbl{image_i} '_image.nii'];
        try 
        img_roi=saveroi(img_mat, lbl_roi);
        save_as_image(img_roi,lbl_roi);
        end
        end
        end
    end
    
    %%
else
    reg_tstat=sqrt(results.rsqr)*sqrt(df/(1-results.rsqr));
    results.tstat=reg_tstat;
    try p_value=tcdf(reg_tstat,df,'upper'); results.pvalue=p_value; end
    
    for m = 1:size(results,2)
        figure('pos',[10 10 1200 900]);
        subplot(2,1,1)
        plot(fitstats(m).y,'b');
        hold on; plot(fitstats(m).yhat,'r')
        dataNames = {'Actual data','predicted data'};
        legend(dataNames)
        rf='Pooled data from all subjects';
        title(rf, 'interpreter','none')
        ylabel(sprintf('R2 = %.2f',results(m).rsqr));
        xlabel(['BIC=' num2str(results(m).bic) '  AIC=' num2str(results(m).aic)])
        subplot(2,1,2)
        bar(results(m).b(1:size(results(1,m).ROIs,2)))
        title([modelName ' - parameters, Intercepts for each subject are not shown'], 'interpreter','none')
        lbl=cell(1,size(results.ROIs,2));
        for x_lbl=1:size(results.ROIs,2)
            lbl{1,x_lbl}=[results.ROIs{1,x_lbl} ' at ' results.timings{1,x_lbl}];
            try
                if results.pval(x_lbl,1) <0.05
                    lbl{1,x_lbl}=[lbl{1,x_lbl} '*'];
                end
            end
        end
        set(gca,'XTick',1:size(results.ROIs,2));
        try xlabel_oblique(lbl); end
    end
end

    fprintf('Saving figures \n')
    h=findobj('type','figure'); % find the handles of the opened figures
    folder=fit.dest_dir;  % Desination folder
    for k=1:numel(h)
        filename=sprintf('%d.jpg',k);
        filename2=sprintf('%d.fig',k);
        fileset=name_fitBOLD(fit);
        figure_name=[fileset '_' filename];
        figure_name2=[fileset '_' filename2];
        file=fullfile(folder,figure_name);
        file2=fullfile(folder,figure_name2);
        if ~isempty(h)
            saveas(h(k),file)
            saveas(h(k),file2)
        end
        close
    end

end
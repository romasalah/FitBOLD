function [rescell,sig_results,rescell_Bmodel]=orgres_fitBOLD(results,fitstats,fit,sub,pooled)
if pooled==0 || isempty(strfind(fit.dependence,'fused'))
    sig_results=results;
    for sub_id=1:size(results,1)
        if isempty(find(results(sub_id).pval<0.05))
            sig_results(sub_id).aic=0;
            sig_results(sub_id).bic=0;
            sig_results(sub_id).train_yhat=mean(results(sub_id).y,1)+zeros(size(results(sub_id).y,1),1);
            sig_results(sub_id).betas=zeros(size(results(sub_id).betas,1),1);
            sig_results(sub_id).ROIs='No significant terms';
            sig_results(sub_id).timings=' ';
        end
        if isempty(strfind(fit.dependence,'step')) && isempty(strfind(fit.dependence,'OLS'))
            if fitstats{sub_id,1}.model_pval>0.05
                sig_results(sub_id).model_in=-1;
            else
                sig_results(sub_id).model_in=1;
            end
        end
        rescell{sub_id,1}=sub(sub_id);
        rescell{sub_id,2}=results(sub_id).rsqr;
        rescell{sub_id,4}=results(sub_id).bic;
        rescell{sub_id,5}=results(sub_id).betas;
        rescell{sub_id,6}=results(sub_id).ROIs;
        rescell{sub_id,7}=results(sub_id).pval;
        rescell{sub_id,8}=results(sub_id).timings;
        rescell{sub_id,9}=fitstats{sub_id,1};
        rescell{sub_id,11}=results(sub_id,1).nlag;
        if ~isempty(strfind(fit.dependence,'fused'))
            rescell{sub_id,12}=fitstats{sub_id,1}.lambda1;
            rescell{sub_id,13}=fitstats{sub_id,1}.lambda2;
            rescell{sub_id,14}=fitstats{sub_id,1}.lambda3;
        end
        if fit.AR~=0
            rescell{sub_id,15}=fitstats{sub_id,1}.lagged_BOLDr2_fused;
            rescell{sub_id,16}=fitstats{sub_id,1}.lagterms_r2;
            rescell{sub_id,17}=fitstats{sub_id,1}.lagged_BOLD_bic;
            rescell{sub_id,18}=fitstats{sub_id,1}.lagterms_bic;
            rescell{sub_id,10}=fitstats{sub_id,1}.rbar;
        end
        if isempty(strfind(fit.dependence,'step'))
            if ~isempty({rescell{sub_id,6}})
                for comb_i=1:size(rescell{sub_id,6},2)
                    rescell{sub_id,19}{1,comb_i}=[rescell{sub_id,6}{1,comb_i} ' at ' rescell{sub_id,8}{1,comb_i}];
                end
            end
        else
            rescell{sub_id,19}=results(sub_id).ROIs;
        end
    end
    rescell_Bmodel=0;
else 
    rescell=cell(1,20);
    rescell_Bmodel=cell(1,20);
    for m = 1:size(results,2)
        
        rescell{m,1}=m;               rescell_Bmodel{m,1}=m;
        rescell{m,2}=results(m).rsqr;   rescell_Bmodel{m,2}=fitstats{m,1}.rsqr;
        rescell{m,3}=results(m).aic; rescell_Bmodel{m,3}=fitstats{m,1}.aic ;
        rescell{m,4}=results(m).bic;  rescell_Bmodel{m,4}=fitstats{m,1}.bic;
        rescell{m,5}=results(m).betas; rescell_Bmodel{m,5}=fitstats{m,1}.beta;
        rescell{m,6}=results(m).ROIs;  rescell_Bmodel{m,6}=results(model_b).ROIs; %stable in the best model
        rescell{m,7}=results(m).pval'; rescell_Bmodel{m,7}=fitstats{m,1}.pval'; %Transposing this pval vector makes a confusion
        rescell{m,8}=results(m).timings;rescell_Bmodel{m,8}=results(model_b).timings; %stable in the best model
        rescell{m,9}=results(m).prediction_r2;    rescell_Bmodel{m,9}=fitstats{m,1}.prediction_r2;
        rescell{m,10}=results(m).prediction_pval; rescell_Bmodel{m,10}=fitstats{m,1}.prediction_pval;
        
        rescell{m,12}=results(m).lambda1; %stable in the best model
        rescell{m,13}=results(m).lambda2;%stable in the best model
        rescell{m,14}=results(m).lambda3;%stable in the best model
        rescell{m,15}=results(m).lambda_BIC;%not done in the best model
        rescell{m,17}=results(m).model_in;%it's only one model
        if AR~=0
            rescell{m,11}=results(m).nlag;            rescell_Bmodel{m,11}=fitstats{m,1}.nlag;
            rescell{m,21}=results(m).prediction_r2_lagcorct;   rescell_Bmodel{m,21}=fitstats{m,1}.prediction_r2_lagcorct;
            rescell{m,22}=results(m).prediction_pval_lagcorct; rescell_Bmodel{m,22}=fitstats{m,1}.prediction_pval_lagcorct;
        end
        if ~isempty({rescell{m,6}}) && ~isempty({rescell_Bmodel{m,6}})
            for comb_i=1:size(rescell{m,6},2)
                rescell{m,19}{1,comb_i}=[rescell{m,6}{1,comb_i} ' at ' rescell{m,8}{1,comb_i}];
            end
            for comb_i2=1:size(rescell_Bmodel{m,6},2)
                rescell_Bmodel{m,19}{1,comb_i2}=[rescell_Bmodel{m,6}{1,comb_i2} ' at ' rescell_Bmodel{m,8}{1,comb_i2}];
            end
        end
        
    end
    sig_results=0;
end
end
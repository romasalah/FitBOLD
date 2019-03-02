function fit_settings=Batch_fitBOLD(fit)
for glm_i=1:size(fit.GLMs,2)
    for pooled_i= 1:size(fit.pooled,2)
        for models_i=1:size(fit.models,2)
            for soc_typ=1:size(fit.social,2)
                fit.condition_before=fit.social{1,soc_typ};
                for social_ROIs_i=1:size(fit.social_ROIs,2)
                    for shifts_i=1:size(fit.shifts,2)
                        for combine_i=1:size(fit.combine,2)
                            for zscore_i= 1:size(fit.zscored,2)
                                for use_subROIs_i=1:size(fit.use_subROIs,2)
                                curfit=fit; curfit.use_subROIs=fit.use_subROIs(1,use_subROIs_i);
                                curfit.GLMs= fit.GLMs{1,glm_i}; curfit.models=fit.models{1,models_i};
                                curfit.condition_before=fit.condition_before;  curfit.shifts=fit.shifts{shifts_i};
                                curfit.combine=fit.combine(1,combine_i); curfit.pooled=fit.pooled(1,pooled_i);
                                curfit.zscored=fit.zscored(zscore_i); curfit.social_ROIs=fit.social_ROIs(social_ROIs_i);
                                main_fitBOLD(curfit);
                                close all %keep this to avoid hundreds of figures on your screen
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
fit_settings=fit;
end

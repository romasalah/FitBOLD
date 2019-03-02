function comparemodel_fitBOLD(folder,type)
if folder==0
    folder=pwd;
end
theories={'IS','RL'};
    conds={'active','passive','n-soc'};
    cd(folder);
    for curt_i=1:size(theories,2)
        curth=theories{curt_i};
    query= ['*' curth '*' type '*results.mat'];

    curt(curt_i).files= dir(query);
    for load_i=1: size( curt(curt_i).files,1)
        try  curt(curt_i).results(load_i)=load( curt(curt_i).files(load_i).name); catch  curt(curt_i).failed_loads; end
    end
    names=erase_fitBOLD( curt(curt_i).files,query,'',conds);   
        for conds_i=1:size( curt(curt_i).results,2)
            curt(curt_i).res(conds_i)=   curt(curt_i).results(conds_i).results;
            curt(curt_i).onestr{conds_i,1}=curt(curt_i).res(conds_i).summaryci.adjrsqr;
            curt(curt_i).onestr{conds_i,2}=curt(curt_i).res(conds_i).summaryci.RMSE;
            curt(curt_i).onestr{conds_i,3}=curt(curt_i).res(conds_i).summaryci.AIC;
            curt(curt_i).onestr{conds_i,4}=curt(curt_i).res(conds_i).summaryci.rsqr;
            curt(curt_i).onestr{conds_i,5}=curt(curt_i).res(conds_i).summaryci.BIC;
            curt(curt_i).onestr{conds_i,6}=curt(curt_i).res(conds_i).summaryci.coefnum;
            curt(curt_i).fitstats{conds_i,1}=curt(curt_i).res(conds_i).summaryci.fstat;
            curt(curt_i).fitstats{conds_i,2}=curt(curt_i).res(conds_i).summaryci.pval(1);
        end
   lbl={names.name};
    for conds_i=1:size(conds,2)
        curt(curt_i).conds{2,conds_i}=find(~cellfun(@isempty,strfind(lbl,conds{1,conds_i})));
    end
    curt(curt_i).lbl_ord=[curt(curt_i).conds{2,1:end}];
    
    curt(curt_i).catsocial=curt(curt_i).onestr;
    for res_i= 1:size(curt(curt_i).onestr,2)
                curt(curt_i).catsocial{curt(curt_i).lbl_ord(1),res_i}= [curt(curt_i).catsocial{curt(curt_i).lbl_ord(1),res_i} curt(curt_i).catsocial{curt(curt_i).lbl_ord(2),res_i} curt(curt_i).catsocial{curt(curt_i).lbl_ord(3),res_i}];
                curt(curt_i).catsocial{curt(curt_i).lbl_ord(2),res_i}=[];
                curt(curt_i).catsocial{curt(curt_i).lbl_ord(3),res_i}=[]; 
    end
    curt(curt_i).catsocial2={};
    for o=1:size(curt(curt_i).lbl_ord,2)
        curt(curt_i).catsocial2=vertcat(curt(curt_i).catsocial2,{curt(curt_i).catsocial{curt(curt_i).lbl_ord(o),:}});
    end
    curt(curt_i).catsocial=curt(curt_i).catsocial2;
 
        curt(curt_i).nonempt_i= find(~cellfun(@isempty,{curt(curt_i).catsocial{1:end,1}}));
    curt(curt_i).catplot_array=curt(curt_i).catsocial([curt(curt_i).nonempt_i],:);
     end  
    
    %final array and plotting
    colors={'r','y','g'};
     catplot_arraymu=cell(1,size(curt(1).catplot_array,2)); 
     catplot_arrayhi=cell(1,size(curt(1).catplot_array,2));  
     catplot_arraylo=cell(1,size(curt(1).catplot_array,2)); 
    for model_i=1:size(theories,2)
        for vi=1:size(curt(model_i).catplot_array,2)
            catplot_arraymu{vi}=horzcat( catplot_arraymu{vi},curt(model_i).catplot_array{1,vi}(1,:)');
            catplot_arraylo{vi}=horzcat( catplot_arraylo{vi},curt(model_i).catplot_array{1,vi}(2,:)');
            catplot_arrayhi{vi}=horzcat( catplot_arrayhi{vi} ,curt(model_i).catplot_array{1,vi}(3,:)');
        end
    end
        ctlbl='';
   
    titles={'Adjusted Rsquared','RMSE','AIC','Ordinary Rsquared','BIC','Number of ROIs'};
        for ct_i=1:size(titles,2)
         titles{1,ct_i}=[  titles{1,ct_i} ctlbl];
        end
     
        for i=1:6
            figure;
            bar(catplot_arraymu{i})
        end
        

    subperplot=3;
    for fig_i=1:2
        figure('Name',['Model fit comparsion' num2str(fig_i) ctlbl], 'Position',[20,20, 800,1200]);
        ni=(fig_i-1)*subperplot;
        for p_i=1:subperplot
            subp=subplot(3,1,p_i);
            subp.XLim=[0 size(conds,2)+1];
            for conds_i=1:size(conds,2)
                mu=zeros(size(conds,2),size(theories,2));
                mu(conds_i,:)=catplot_arraymu{ni+p_i}(conds_i,:);
                loci=catplot_arraylo{ni+p_i}(conds_i,:);
                hici=catplot_arrayhi{ni+p_i}(conds_i,:);
                n4{conds_i}=bar(mu,colors{1,conds_i});  hold on;
                
%% write the value on top of bars
                for i=1:length(n4{conds_i})
                    xpos(i,:)=n4{conds_i}(i).XData; xpos2bars=xpos(i,conds_i);
                    ypos(i,:)=n4{conds_i}(i).YData; ypos2=ypos(i,conds_i);
                    ydat(i,:)=n4{conds_i}(i).YData; ydat2=ydat(i,conds_i);
                    
                    switch i; case 1; xpos2barstt=xpos2bars-0.15; case 2; xpos2barstt=xpos2bars+0.15; end
                    if ydat2>10; format='%.d'; else; format='%.3f';  end
                   if ydat2~=0
                    text(xpos2barstt,ypos2,num2str(ydat2,format),'VerticalAlignment','bottom','horizontalalign','center');
                    if i==2
                        text(xpos2barstt,ypos2+ypos2*0.2,conds{1,conds_i},'VerticalAlignment','bottom','horizontalalign','center');
                    end
                    end
                xpossubp(1,conds_i)=xpos2bars;
                end
                clear xpos2bars
                
            end
  %% extend the Axes          
            allvals=[catplot_arrayhi{1:end,ni+p_i} catplot_arraylo{1:end,ni+p_i}];
                if mean(allvals,2)>0; maxval=max(allvals(:)); else maxval=min(allvals(:)); end
                 ylims=sort([ 0  maxval+30*maxval/100]);
                 subp.YLim=ylims;
           
    %% reposition and label the axes        
            title(titles{1,ni+p_i})
            lo=[];
            xpossubp2bars=[xpossubp-0.15; xpossubp+0.15];lo=xpossubp2bars;
            xpossubp2bars=[lo(:,1); lo(:,2); lo(:,3)]';
            subp.XTick=xpossubp2bars;
            lblobq={};
            roislbl=theories;
            for n=1:size(conds,2)
            lblobq=horzcat(lblobq,roislbl);
            end
            xlabel_oblique(lblobq,70);
            
    %% make error bars        
            for conds_i=1:size(conds,2)
                mu=zeros(size(conds,2),size(theories,2));
                mu(conds_i,:)=catplot_arraymu{ni+p_i}(conds_i,:);
                loci=catplot_arraylo{ni+p_i}(conds_i,:);
                hici=catplot_arrayhi{ni+p_i}(conds_i,:);
                errordist(conds_i,:)=(hici-loci)./2;
                for nbars=1:size(theories,2)
                    allbarloc(nbars,:)=n4{conds_i}(nbars).XData+n4{conds_i}(nbars).XOffset;
                    mybarloc(nbars)=allbarloc(nbars,conds_i);
                    errorbar(mybarloc(nbars),mu(conds_i,nbars),errordist(nbars))
                    clear allbarloc mybarloc
                end
            end
        end
         end
end
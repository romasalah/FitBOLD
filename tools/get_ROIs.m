function ROIs=get_ROIs(drawF)
linpos=orthviews_fitBOLD(0,'FFX_*RL','linear','pos',drawF);
linneg=orthviews_fitBOLD(0,'FFX_*RL','linear','neg',drawF);
 quadpos=orthviews_fitBOLD(0,'FFX_*RL','quad','pos',drawF);
 quadneg=orthviews_fitBOLD(0,'FFX_*RL','quad','neg',drawF);
 for i2=1:9
gg{1,i2}=linpos(i2).model;
gg{2,i2}=linpos(i2).ROIs;
gg{3,i2}=linneg(i2).ROIs;
gg{4,i2}=quadpos(i2).ROIs;
gg{5,i2}=quadneg(i2).ROIs;
 end
gg2=gg';
gg3=cell2table([{gg2{:,2}}' {gg2{:,3}}' {gg2{:,4}}' {gg2{:,5}}']);
gg3.Properties.VariableNames={'Linear_Positive','Linear_Negative', 'Quadratic_Postive','Quadratic_Negative'};
gg3.Properties.RowNames={gg2{1:end,1}};
ROIs=gg3;
save('SoDec_all_ROIs.mat','ROIs');
end
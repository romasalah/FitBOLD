function [termsquad,namesquad,contterms,names,tgt]=Polynom_fitBOLD(design,power)
tgt=design{:,end};
contterms=design{:,1:end-2};
roissize=size(contterms,2);
names={design.Properties.VariableNames{1,1:roissize}};
%% Create polynomials
termsquad=[];namesquad={};
for fill_poly=1:power
    for power_i=1:roissize
        termstmp(:,power_i)=contterms(:,power_i).^fill_poly;
        if fill_poly > 1
            namestmp{1,power_i}=[names{1,power_i} '^' num2str(fill_poly)];
        else
            namestmp{1,power_i}=names{1,power_i};
        end
    end
    termsquad=[termsquad termstmp];
    namesquad=horzcat(namesquad,namestmp);
    clear termstmp namestmp
end
end
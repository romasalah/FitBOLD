function savefig_fitBOLD(name,destination)
fprintf('Saving figures \n')
h=findobj('type','figure'); % find the handles of the opened figures
for k=1:numel(h)
    filename=sprintf('%d.jpg',k);
    filename2=sprintf('%d.fig',k);
    figure_name=[name '_' filename];
    figure_name2=[name '_' filename2];
    file=fullfile(destination,figure_name);
    file2=fullfile(destination,figure_name2);
    if ~isempty(h)
        saveas(h(k),file)
        saveas(h(k),file2)
    end
    close
end
end
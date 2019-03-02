function filter=BOLD_filt(n,c,ar)
max=size(n,1);
order=round(max/3)-2;
if order >50; order=50;end
if exist('c','var') && strcmp(c,'LSD')
    filter=designfilt('Highpassfir',...
  'StopbandFrequency',0.03,'PassbandFrequency',0.035','FilterOrder', order,...
   'DesignMethod','equiripple');
elseif exist('ar','var') && ar>1
filter=designfilt('Bandpassfir',...
  'Cutofffrequency1',0.03,'Cutofffrequency2',1/ar','FilterOrder', order...
   );   
else
filter=designfilt('Bandpassfir',...
  'Cutofffrequency1',0.03,'Cutofffrequency2',0.99','FilterOrder', order...
   );
end
end
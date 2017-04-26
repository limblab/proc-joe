function [] = saveFigure( figHandle, targetDirectory,filename )
% saves figHandle to targetDirectory/... as a PDF, EPS, FIG, and PNG

% check if targetDirectory ends with filesep
if(targetDirectory(end) ~= filesep)
    targetDirectory = [targetDirectory filesep];
end
% check if all folders exist
if exist(strcat(targetDirectory,'Raw_Figures'),'file')~=7
    mkdir(strcat(targetDirectory,'Raw_Figures'))
end
if exist(strcat(targetDirectory,['Raw_Figures' filesep 'PDF']),'file')~=7
    mkdir(strcat(targetDirectory,['Raw_Figures' filesep 'PDF']))
end
if exist(strcat(targetDirectory,['Raw_Figures' filesep 'FIG']),'file')~=7
    mkdir(strcat(targetDirectory,['Raw_Figures' filesep 'FIG']))
end
if exist(strcat(targetDirectory,['Raw_Figures' filesep 'EPS']),'file')~=7
    mkdir(strcat(targetDirectory,['Raw_Figures' filesep 'EPS']))
end
if exist(strcat(targetDirectory,['Raw_Figures' filesep 'PNG']),'file')~=7
    mkdir(strcat(targetDirectory,['Raw_Figures' filesep 'PNG']))
end

% print all file types
filename(filename==' ')='_';%replace spaces in name for saving
print('-dpdf',figHandle,strcat(targetDirectory,['Raw_Figures' filesep 'PDF' filesep],filename,'.pdf'))
print('-deps',figHandle,strcat(targetDirectory,['Raw_Figures' filesep 'EPS' filesep],filename,'.eps'))
print('-dpng',figHandle,strcat(targetDirectory,['Raw_Figures' filesep 'PNG' filesep],filename,'.png'))
saveas(figHandle,strcat(targetDirectory,['Raw_Figures' filesep 'FIG' filesep],filename,'.fig'),'fig')

end


%% 作圖功能
clc;
clear;
close all;
foldername = 'Concentration_Test';
mkdir(foldername);
%eval(['cd ' foldername]);
path=pwd;
Folderpath = [path '\' foldername];


RawLocdata=csvread('Z_Location.txt'); %1=Elenum 2~9=Node Location
LocElem=RawLocdata(:,1);
Location=RawLocdata(:,2:9);
hold off;
c=[0.5 0.5 0.5];


for Iteration=1:112
    Scale=csvread(['D:\Project\Cell_Differentiation_ANSYS\Cell_differentiation_renew_units-m,kg,s\Diffusion\Concentration_' num2str(Iteration) '.inp']);
    for inn=1:length(LocElem)
        patch(Location(inn,1:2:7),Location(inn,2:2:8),c,'FaceAlpha',Scale(inn));
        hold on;
    end
    title(['Iteration' num2str(Iteration) '  ' 'Day ' num2str(Iteration)]);
    saveas(gca,[Folderpath '\Iteration_Result_' num2str(Iteration) '.jpg']);
    close all;
end
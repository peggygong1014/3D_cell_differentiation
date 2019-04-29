function f=plotFunction(DiffType,E_10,Iteration,Folderpath,E_10_Option,type,SPNum,Elenum,resorpElem)
%Type = 1 正常作圖
%Type = 2 不考慮濃度作圖
%Type = 3 Stimulus 作圖
%% 作圖功能
close all;
RawLocdata=csvread('Z_Location.txt'); %1=Elenum 2~9=Node Location
LocElem=RawLocdata(:,1);
Location=RawLocdata(:,2:9);
hold off;
for inn=1:length(LocElem)
    %Resorption
    if sum(resorpElem==Elenum(inn))>0
        DiffType(inn)=0;
    end
    
    c=autoColor(DiffType(inn),type);
    [Perm,P,NONE]=AutoMP_Select(DiffType(int16(inn)));
    Scale = floor(E_10(inn,end-1+E_10_Option)/NONE/0.2);
    Scale = autoScale(Scale,type);
    
    %Resorption
    if DiffType(int16(inn))==0
        Scale=1;
    end
    
    patch(Location(inn,1:2:7),Location(inn,2:2:8),c,'FaceAlpha',Scale,'LineWidth',0.001);
    hold on;
end
title(['Iteration' num2str(Iteration) '  ' 'Day ' num2str(Iteration)]);
axis image;
print(gcf, '-dpng', '-r512',[Folderpath '\Iteration_Result_' SPNum '_' num2str(Iteration,'%04i') '.png']);
close all;
end


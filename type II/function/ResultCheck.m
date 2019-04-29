for i=1:length(OutHistoryType(1,:))-1
    %Type 需逐漸變小
    sum(OutHistoryType(:,i+1)>OutHistoryType(:,i))
end

for i=1:length(OutHistoryType(1,:))-1
    %E需逐漸變大
    %變化後Young's  變化前Young's
    sum(OutE(:,i+1)<OutE(:,i))
end

i=1;
checkNum=2;
checkNum=checkNum+1;
Enum=find(OutE(:,i+1)<OutE(:,i));

ResultC=E_10(Enum(checkNum),boolean([0 E_10(Enum(checkNum),2:end-1)>0 0]));

for Test=1:length(ResultC)-1
    if(ResultC(Test+1)<ResultC(Test))
        disp(Test)
    end
end

for Test=1:length(ResultC)
    disp(mean(ResultC(1:Test)))
end

OutHistoryType(checkNum,:)

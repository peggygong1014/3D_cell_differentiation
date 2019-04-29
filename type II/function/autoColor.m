function [ c ] = autoColor( DiffType, type )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

switch DiffType
    %Resorption
    case 0
        c=[255 255 255];
    case 1 %matture bone(Cancellous bone) 0.266>=S
        %             c=[255 80 0];
        c=[0 255 0];
        if (type==3)
            c=[78 255 51];
        end
    case 2 %immature bone 1>=S>0.266
        %             c=[255 133 133];
        c=[200 255 200];
        if (type==3)
            c=[0 163 16];
        end
    case 3 %cartilage 3>=S>1
        %             c=[255 180 0];
        c=[90 10 255];
        if (type==3)
            c=[255 162 0];
        end
    case 4 %fibrous connective tissue S>3
        %             c=[0 178 255];
        c=[255 0 255];
        if (type==3)
            c=[255 255 0];
        end
end
c=c./255;
end


function [ Scale ] = autoScale( Scale,type )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    switch Scale
        case 0
%             Scale = 0.1;
            Scale = 1;
            if (type ==2||type==3)
                Scale = 1;
            end
        case 1
%             Scale = 0.2;
            Scale = 1;
            if (type ==2||type==3)
                Scale = 1;
            end
        case 2
%             Scale = 0.32;
            Scale = 1;
            if (type ==2||type==3)
                Scale = 1;
            end
        case 3
%             Scale = 0.5;
            Scale = 1;
            if (type ==2||type==3)
                Scale = 1;
            end
        case 4
            Scale = 1;
        case 5
            Scale = 1;
    end
end


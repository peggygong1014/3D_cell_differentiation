function png2movie( path,crticname )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
file=dir([path '\*' crticname '*']);
l=length(file);
writerObj = VideoWriter([path '\Result_animation.avi'],'Uncompressed AVI');
writerObj.FrameRate = 2;
open(writerObj);
for i=1:l
    A = imread([path '\' file(i).name]);
    writeVideo(writerObj,A)
end
close(writerObj)
end


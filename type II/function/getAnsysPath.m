function ansysPath = getAnsysPath( version )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
ansysPath=[getenv(['ANSYS' version '_DIR']) '\bin\winx64\ANSYS' version '.exe'];
end


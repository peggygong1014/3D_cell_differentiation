function preParameter
global E_granulation E_fibrous E_cartilage E_marrow E_immature E_mature E_cortical ansysPath
global Perm_granulation Perm_fibrous Perm_cartilage Perm_marrow Perm_immature Perm_mature Perm_cortical
global P_granulation P_mature P_cortical P_fibrous P_cartilage P_marrow P_immature P_Imp E_Imp APPForce T

ansys_version='170';
APPForce=10; % um
% Material Properties
E_Imp=113000000000;
E_granulation = 1000000;
E_fibrous = 2000000;
E_cartilage = 10000000;
E_marrow = 6000000000;
E_immature = 1000000000;
E_mature = 6000000000;
E_cortical = 20000000000;
P_Imp = 0.3;
Perm_granulation = 1E-14;
Perm_fibrous = 1E-14;
Perm_cartilage = 5E-15;
Perm_marrow = 3.7E-13;
Perm_immature = 1E-13;
Perm_mature = 3.7E-13;
Perm_cortical = 1E-17;
P_granulation = 0.17;
P_fibrous = 0.17;
P_cartilage = 0.17;
P_marrow = 0.3;
P_immature = 0.3;
P_mature = 0.3;
P_cortical = 0.3;
ansysPath=getAnsysPath(ansys_version);
T=1;

fid=fopen('preParameter.inp','w');
fprintf(fid,'APPForce=%E\n',APPForce);
fprintf(fid,'E_Imp=%E\n',E_Imp);
fprintf(fid,'E_granulation=%E\n',E_granulation);
fprintf(fid,'E_fibrous=%E\n',E_fibrous);
fprintf(fid,'E_cartilage=%E\n',E_cartilage);
fprintf(fid,'E_marrow=%E\n',E_marrow);
fprintf(fid,'E_immature=%E\n',E_immature);
fprintf(fid,'E_mature=%E\n',E_mature);
fprintf(fid,'E_cortical=%E\n',E_cortical);
fprintf(fid,'P_Imp=%E\n',P_Imp);
fprintf(fid,'Perm_granulation=%E\n',Perm_granulation);
fprintf(fid,'Perm_fibrous=%E\n',Perm_fibrous);
fprintf(fid,'Perm_cartilage=%E\n',Perm_cartilage);
fprintf(fid,'Perm_marrow=%E\n',Perm_marrow);
fprintf(fid,'Perm_immature=%E\n',Perm_immature);
fprintf(fid,'Perm_mature=%E\n',Perm_mature);
fprintf(fid,'Perm_cortical=%E\n',Perm_cortical);
fprintf(fid,'P_granulation=%E\n',P_granulation);
fprintf(fid,'P_fibrous=%E\n',P_fibrous);
fprintf(fid,'P_cartilage=%E\n',P_cartilage);
fprintf(fid,'P_marrow=%E\n',P_marrow);
fprintf(fid,'P_immature=%E\n',P_immature);
fprintf(fid,'P_mature=%E\n',P_mature);
fprintf(fid,'P_cortical=%E\n',P_cortical);
fprintf(fid,'T=%E\n',T);
fclose(fid);
copyfile('preParameter.inp',[pwd '\Model\']);



end
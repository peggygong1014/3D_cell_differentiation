clc
clear
global E_granulation E_fibrous E_cartilage E_marrow E_immature E_mature E_cortical
global Perm_granulation Perm_fibrous Perm_cartilage Perm_marrow Perm_immature Perm_mature Perm_cortical
global P_granulation P_mature P_cortical P_fibrous P_cartilage P_marrow P_immature

%% 參數設置
format longEng
teststep=28;
OutputMPMAX=300;
E_10_Option=1; %10步平均開關 1=使用 0=不使用
oneDirectionPathOption=0; %分化單向開關,不可逆1=使用 0=不使用
digitalOption=0;
a=0.0375;
b=3/1E6; %m/s to um/s (3um/s)xc.4u
APPForce=10; % um

% Material Properties


E_Imp=113000000000;
E_granulation = 200000;
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

% FluidHist=csvread(['D:\Project\Cell_Differentiation_ANSYS\Cell_differentiation_renew_units-m,kg,s\2016-9-22 16.14_500N\FluidHist.csv']);
% StrainHist=csvread(['D:\Project\Cell_Differentiation_ANSYS\Cell_differentiation_renew_units-m,kg,s\2016-9-22 16.14_500N\StrainHist.csv']);
foldername = [clock];
foldername = ['Strain_Fluid_Test\' num2str(foldername(1)) '-' num2str(foldername(2)) '-' num2str(foldername(3)) ' ' num2str(foldername(4)) '.' num2str(foldername(5)) ];

mkdir(foldername);
path=pwd;
Folderpath = [path '\' foldername];

copyfile('2.Solve.inp',Folderpath);
copyfile('Z_Model_Build.inp',Folderpath);
copyfile('Z_Strain_Fluid_Test.m',Folderpath);




a=0.0375;
b=3/1E6; %m/s to um/s (3um/s)

Iteration=1;
OperatingResult=system(['SET KMP_STACKSIZE=2048k & "D:\ANSYS Inc\v150\ansys\bin\winx64\ANSYS150.exe" -b -j "MATLAB_Differentiation" -i Z_Model_Build.inp -o Output MATLAB_PART'])
OperatingResult=system(['SET KMP_STACKSIZE=2048k & "D:\ANSYS Inc\v150\ansys\bin\winx64\ANSYS150.exe" -b -j "MATLAB_Differentiation" -i 1.Presetting.inp -o Output MATLAB_PART']);
OperatingResult=system(['SET KMP_STACKSIZE=2048k & "D:\ANSYS Inc\v150\ansys\bin\winx64\ANSYS150.exe" -b -j "MATLAB_Differentiation" -i 2.Solve.inp -o Output MATLAB_PART'])

Concentration=csvread([path '\Diffusion\Concentration_' num2str(Iteration) '.inp']);
RawData = csvread('Z_PredictMat.txt');
Elenum=RawData(:,1);
Prin1 = RawData(:,2);
Prin2 = RawData(:,3);
Prin3 = RawData(:,4);
FlowX = RawData(:,5);
FlowY = RawData(:,6);
Flow = RawData(:,7);
EPX = RawData(:,8);
EPY = RawData(:,9);
EPZ = RawData(:,10);
EPXY = RawData(:,11);
EPYZ = RawData(:,12);
EPXZ = RawData(:,13);

%Differentiation公式計算
Octahedral_strain = sqrt((Prin1-Prin2).^2+(Prin2-Prin3).^2+(Prin3-Prin1).^2)./2; %Modeling Mechanical Signals on the Surface ofmCT and CAD Based Rapid Prototype ScaffoldModels to Predict (Early Stage) TissueDevelopment
% Octahedral_strain = sqrt((EPX-EPY).^2+(EPY-EPZ).^2+(EPZ-EPX).^2+6.*((EPXY./2).^2+(EPYZ./2).^2+(EPXZ./2).^2)).*(2/3); %http://www.colorado.edu/MCEN/MCEN5023/chap_04.pdf
% Octahedral_strain = sqrt((Prin1-Prin2).^2+(Prin2-Prin3).^2+(Prin3-Prin1).^2).*(2/3); %https://books.google.com.tw/books?id=PDF7INRk1tsC&pg=PA272&lpg=PA272&dq=Octahedral+shear+strain&source=bl&ots=UVKAMyIJEO&sig=mpQ1zDP9Mfav4lFSWFceG0YGzN0&hl=zh-TW&sa=X&ved=0ahUKEwiqzs22xtnPAhUMuo8KHQv6D3wQ6AEIfjAO#v=onepage&q=Octahedral%20shear%20strain&f=false
% Octahedral_strain = abs(4.*Prin1.*Prin2.*Prin3).^(1/3); %Distortion stress, distortion strain and their physical concept
% Octahedral_strain = sqrt(Prin1.^2-(Prin1.*Prin2)+(Prin2).^2).*sqrt(2/3); % A prediction of cell differentiation and proliferation within a collagen–glycosaminoglycan scaffold subjected to mechanical strain and perfusive fluid flow
% Octahedral_strain = sqrt((Prin1-Prin2).^2+(Prin1-Prin3).^2+(Prin2-Prin3).^2)./sqrt(2); %A numerical model of the fracture healing process that describes tissue development and revascularisation

% Octahedral_strain = abs(Prin1-Prin2)./2;
% Octahedral_strain = sqrt(  ((EPX-EPY).^2+(EPY-EPZ).^2+(EPZ-EPX).^2)./9  +  (EPXY.^2+EPYZ.^2+EPXZ.^2)./6  ); %Website
% Octahedral_strain = sqrt((EPX-EPY).^2+(EPY-EPZ).^2+(EPZ-EPX).^2).*(3/2);

Fluid_velocity = sqrt(abs((FlowX.^2 + FlowY.^2)));
E_10_Option=1;

for ip=1:2
        if (ip== 1)
            SPNum='Strain';
            a=0.0375;
            b=999;
        else
            SPNum='Fluid';
            a=999;
            b=3/1E6; %m/s to um/s (3um/s)
        end
        % TestS=(StrainHist./a + FluidHist./b);
        TestS=(Octahedral_strain./a + Fluid_velocity./b);
        % Cell_type: 0=granulation tissue 1=mature bone 0.266>=S 2=immature bone 1>=S>0.266 3=cartilage 3>=S>1 4=fibrous connective tissue S>3
        TestDiffType = zeros(size(TestS));
        TestDiffType(TestS>3) = 4;
        TestDiffType(TestS>1 & TestS<=3) = 3;
        TestDiffType(TestS>0.266 & TestS<=1) = 2;
        TestDiffType(TestS<=0.266) = 1;
        E_10Test=ones(length(Octahedral_strain),1).*2000000;
        plotFunction(TestDiffType(:,1),E_10Test,1,Folderpath,E_10_Option,3,SPNum);
    end
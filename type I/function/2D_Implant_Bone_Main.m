clc;
clear;

global E_granulation E_fibrous E_cartilage E_marrow E_immature E_mature E_cortical
global Perm_granulation Perm_fibrous Perm_cartilage Perm_marrow Perm_immature Perm_mature Perm_cortical
global P_granulation P_mature P_cortical P_fibrous P_cartilage P_marrow P_immature

%% �ѼƳ]�m
format longEng
teststep=28;
OutputMPMAX=300;
E_10_Option=1; %10�B�����}�� 1=�ϥ� 0=���ϥ�
oneDirectionPathOption=0; %���Ƴ�V�}��,���i�f1=�ϥ� 0=���ϥ�
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


%% �Ыظ�Ƨ�
foldername = [clock];
foldername = [num2str(foldername(1)) '-' num2str(foldername(2)) '-' num2str(foldername(3)) ' ' num2str(foldername(4)) '.' num2str(foldername(5)) '_' num2str(APPForce) 'N'];
mkdir(foldername);
%eval(['cd ' foldername]);
path=pwd;
Folderpath = [path '\' foldername];

%%
delete 'Term.db';
copyfile('Term-BAK.db','Term.db');
copyfile('A_Main.m',Folderpath);


OutHistoryType=[];
OutE=[];
StrainHist=[];
FluidHist=[];
E_10=[]; %1=ElemNum 2~11=E 12=AverageResult
P_10=[]; %1=ElemNum 2~11=E 12=AverageResult
Perm_10=[]; %1=ElemNum 2~11=E 12=AverageResult

system('SET'); %���F�� MATLAB �i�H�] ANSYS
clc;
disp('Creating Folder Done ')
disp('ParameterSetting Done ')

for Iteration=1:teststep
    
    disp(['----------------Iteration : ' num2str(Iteration) ' Start----------------']);
    %% Pre-Parameter
    fid=fopen(['Z_MATLAB_parameter.inp'],'w');
    fprintf(fid,'teststep = %d\n',Iteration);
    fprintf(fid,'APPForce = %d\n',APPForce);
    fclose(fid);
    %% ANSYS Part
    if(Iteration==1)
    %File 1.Presetting
    OperatingResult=system(['SET KMP_STACKSIZE=2048k & "D:\ANSYS Inc\v150\ansys\bin\winx64\ANSYS150.exe" -b -j "MATLAB_Differentiation" -i 1.Presetting.inp -o Output MATLAB_PART']);
    if OperatingResult~=8
        disp(['File 1.Presetting .....Error']);
        break;
    end
    disp(['File 1.Presetting .....Done']);
    end
    
    %File 2.Solve
    OperatingResult=system(['SET KMP_STACKSIZE=2048k & "D:\ANSYS Inc\v150\ansys\bin\winx64\ANSYS150.exe" -b -j "MATLAB_Differentiation" -i 2.Solve.inp -o Output MATLAB_PART']);
    if OperatingResult~=8
        disp(['File 2.Solve .....Error']);
        break;
    end
    disp(['File 2.Solve .....Done']);
    
    
    
    %% ���G�Ʀ�� Ex 1.1 �� 1 & 1.05��1 �Ӵ�֧��ưѼƪ�����
    %�ƾڰѼƤ��t
    Concentration=csvread([path '\Diffusion\Concentration_' num2str(Iteration) '.inp']);
    RawData = csvread('Z_PredictMat.txt');
    Elenum=RawData(:,1);
    Prin1 = RawData(:,2);
    Prin2 = RawData(:,3);
    Prin3 = RawData(:,4);
    FlowX = RawData(:,5);
    FlowY = RawData(:,6);
    Flow = RawData(:,7);
    
    %Differentiation�����p��
    Octahedral_strain = sqrt((Prin1-Prin2).^2+(Prin2-Prin3).^2+(Prin3-Prin1).^2)./2;
%     Octahedral_strain = abs(4.*Prin1.*Prin2.*Prin3).^(1/3);
    Fluid_velocity = sqrt(abs((FlowX.^2 - FlowY.^2)));
    Fluid_velocity2 = Flow;
    S=(Octahedral_strain./a + Fluid_velocity./b);
    % Cell_type: 0=granulation tissue 1=mature bone 0.266>=S 2=immature bone 1>=S>0.266 3=cartilage 3>=S>1 4=fibrous connective tissue S>3
    DiffType = zeros(size(S));
    DiffType(S>3) = 4;
    DiffType(S>1 & S<=3) = 3;
    DiffType(S>0.266 & S<=1) = 2;
    DiffType(S<=0.266) = 1;
    
    if Iteration>1
        DiffAvoidInd=DiffType>HistoryType;
        if oneDirectionPathOption == 1
            %�u���\4.Fibrous��3.Cartilage��2.Immature��1.Mature
            DiffAvoidInd=DiffType>HistoryType;
            DiffType(DiffAvoidInd)=HistoryType(DiffAvoidInd);
        end
    end
    
    Result=[];
    Result1=[];
    Result2=[];
    for in=1:length(DiffType)
        [Perm,P,E_tissue]=AutoMP_Select(DiffType(in));
        Result=[Result ; (1-Concentration(in))*E_granulation+Concentration(in)*E_tissue];
        Result1=[Result1 ; (1-Concentration(in))*P_granulation+Concentration(in)*P];
        Result2=[Result2 ; (1-Concentration(in))*Perm_granulation+Concentration(in)*Perm];
    end
    
    
    
    if length(E_10)~=length(Elenum)
        E_10=zeros(length(Elenum),12);
        P_10=zeros(length(Elenum),12);
        Perm_10=zeros(length(Elenum),12);
    end
    if sum(E_10(:,1))==0
        E_10(:,1)=Elenum;
        P_10(:,1)=Elenum;
        Perm_10(:,1)=Elenum;
    end
    
    if Iteration>1
        resetLoc=find(DiffAvoidInd);
        if(isempty(resetLoc)==0)
            E_10(resetLoc,2:end-1)=0;
        end
    end
    
    
    
    filter=Result;
    filter1=Result1;
    filter2=Result2;
    if(digitalOption==1)
        %-------------------------�Ʀ�ƥ\��---------------------------------
        UniMP=unique(Result); %�ƾھ�z�A�h��Raw data���ƪ��ƭ�
        filter=round(filter/filterIndex)*filterIndex; %Filter Data
        Ori=length(UniMP);
        %----------------------------------------------------------
    end
    
    %�������v�����Y�Ƭ��� (���F�Q�B����)
    E_10(:,2:end-2)=E_10(:,3:end-1); %�C�����v�ƾګe��
    P_10(:,2:end-2)=P_10(:,3:end-1); %�C�����v�ƾګe��
    Perm_10(:,2:end-2)=Perm_10(:,3:end-1); %�C�����v�ƾګe��
    %
    E_10(:,end-1)=filter; %�N�̷s�ƭȩ�b�̫�@��
    P_10(:,end-1)=filter1; %�N�̷s�ƭȩ�b�̫�@��
    Perm_10(:,end-1)=filter2; %�N�̷s�ƭȩ�b�̫�@��
    %
    Eaverage = sum(E_10(:,2:end-1),2) ./ sum(E_10(:,2:end-1)~=0,2); %��ثe�O���X�Ӿ��v�ƾ�
    Paverage = sum(P_10(:,2:end-1),2) ./ sum(P_10(:,2:end-1)~=0,2); %��ثe�O���X�Ӿ��v�ƾ�
    Permaverage = sum(Perm_10(:,2:end-1),2) ./ sum(Perm_10(:,2:end-1)~=0,2); %��ثe�O���X�Ӿ��v�ƾ�
    %
    E_10(:,end)=Eaverage; %�x�}�Gsum(a)�O�C�D�M�Asum(a,2)�O��D�M ���B�D�X�Q�B�����é�b�̫�@���
    P_10(:,end)=Paverage; %�x�}�Gsum(a)�O�C�D�M�Asum(a,2)�O��D�M ���B�D�X�Q�B�����é�b�̫�@���
    Perm_10(:,end)=Permaverage; %�x�}�Gsum(a)�O�C�D�M�Asum(a,2)�O��D�M ���B�D�X�Q�B�����é�b�̫�@���
    
    
    UniMP=E_10(:,end-1+E_10_Option);
    UniP=P_10(:,end-1+E_10_Option);
    UniPerm=Perm_10(:,end-1+E_10_Option);
    
    if(digitalOption==1)
        %-------------------------�Ʀ�ƥ\��---------------------------------
        UniMP=[];
        UniMP=unique(E_10(:,end-1+E_10_Option)); %After filt and sort
        disp([ 'Origin Size = ' num2str(Ori) '; After filter = ' num2str(length(UniMP)) ]);
        %----------------------------------------------------------
    end
    
    ProcessResult={}; %�Ĥ@��O���ưѼơA�᭱�d��O�ϥθӧ��ưѼƪ������s��
    
    switch E_10_Option
        case 0
        processColumn='filter';
        case 1
        processColumn='Eaverage';
    end
    switch digitalOption
        case 0
            for i=1:length(UniMP)
                ProcessResult{i}=[UniMP(i) UniP(i) UniPerm(i) Elenum(i)];
            end
        case 1
            for i=1:length(UniMP)
                MPEnum =Elenum(find(eval(processColumn) == UniMP(i)))';
                ProcessResult{i}=[UniMP(i) MPEnum];
            end
    end
            
            
    %% ��X���ưѼƤΤ����s��
    part=1;
    fid=fopen(['MP_Iteration_' num2str(part) '.inp'],'w');
    fprintf(fid,'resume,Term,db\n');
    fprintf(fid,'/prep7\n');
    fprintf(fid,'allsel\n');
    fprintf(fid,['/INPUT,Z_PreSettingClear,''inp'',,, 0\n']);
    mpCount=0;
    for in=1:length(ProcessResult)
        E=ProcessResult{in}(1); % 1=Young's modulus
        P=ProcessResult{in}(2); % 2=Poisson ratio
        Perm=ProcessResult{in}(3); % 3=Permeability
        MPEnum =ProcessResult{in}(4); %��X�������Young's modulus�������s��
        mpCount=mpCount+1; %��e��J�����ưѼƽs��
        fprintf(fid,'mp,EX,%d,%E\n',mpCount+4,E);
        fprintf(fid,'mp,PRXY,%d,%E\n',mpCount+4,P);
        fprintf(fid,'TB,PM,%d,,,PERM\n',mpCount+4);
        fprintf(fid,'TBDATA,1,%E,%E,%E\n',Perm,Perm,Perm);
        for OutputElemNum = 1:length(MPEnum)
            fprintf(fid,'emodif,%d,mat,%d\n',MPEnum(OutputElemNum),mpCount+4);
        end
        
        if mod(mpCount,OutputMPMAX)==0
            fprintf(fid,'save,Term,db\n');
            fclose(fid);
            if in~=length(ProcessResult)
                part=part+1;
                fid=fopen(['MP_Iteration_' num2str(part) '.inp'],'w');
                fprintf(fid,'resume,Term,db\n');
                fprintf(fid,'/prep7\n');
                fprintf(fid,'allsel\n');
            end
        end
    end
    if mod(mpCount,OutputMPMAX)~=0
        fprintf(fid,'save,Term,db\n');
        fclose(fid);
    end
    disp(['Material Properties Output ' num2str(part) '.....Done']);
    %% Plot Result
    plotFunction(DiffType,E_10,Iteration,Folderpath,E_10_Option);
    disp('Plot .....Done');
    
    
    %% Save Data
    HistoryType=DiffType;
    OutHistoryType=[OutHistoryType HistoryType];
    OutE=[OutE E_10(:,end-(1+E_10_Option))];
    
    StrainHist=[StrainHist Octahedral_strain];
    FluidHist=[FluidHist Fluid_velocity];
    csvwrite('HistoryDiffType.csv',OutHistoryType);
    csvwrite('HistoryEavg.csv',OutE);
    movefile('HistoryDiffType.csv',Folderpath,'f');
    movefile('HistoryEavg.csv',Folderpath,'f');
    
    
    %% Input Material Properties
    for MPpart=1:part
        system(['SET KMP_STACKSIZE=2048k & "D:\ANSYS Inc\v150\ansys\bin\winx64\ANSYS150.exe" -b -j "MATLAB_Differentiation" -i MP_Iteration_' num2str(MPpart) '.inp -o Output MP_input']);
        disp(['Input Material Properties ' num2str(MPpart) '  .....Done']);
    end
    
    disp(['----------------Iteration : ' num2str(Iteration) ' Finish---------------']);
    
end
csvwrite('StrainHist.csv',StrainHist);
csvwrite('FluidHist.csv',FluidHist);
csvwrite('E_10.csv',E_10);
movefile('StrainHist.csv',Folderpath,'f');
movefile('FluidHist.csv',Folderpath,'f');
movefile('E_10.csv',Folderpath,'f');

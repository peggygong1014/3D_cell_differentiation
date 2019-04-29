function Cell_Differentiation
addpath([pwd '\function']);
global E_granulation E_fibrous E_cartilage E_marrow E_immature E_mature E_cortical
global Perm_granulation Perm_fibrous Perm_cartilage Perm_marrow Perm_immature Perm_mature Perm_cortical
global P_granulation P_mature P_cortical P_fibrous P_cartilage P_marrow P_immature P_Imp E_Imp APPForce T
global ansysPath
if isempty(E_granulation)
    preParameter;
end


clc;
%% �ѼƳ]�m
format longEng
teststep=84;
OutputMPMAX=99999; %�@��Material Proeprties.inp�n��J�X�ӧ��ưѼ�
E_10_Option=1; %10�B�����}�� 1=�ϥ� 0=���ϥ�
oneDirectionPathOption=0; %���Ƴ�V�}��,���i�f1=�ϥ� 0=���ϥ�
a=0.0375;
b=3/1E6; %m/s to um/s (3um/s)xc.4u
preParameter;

%% �Ыظ�Ƨ�
foldername = clock;
foldername = [num2str(foldername(1)) '-' num2str(foldername(2)) '-' num2str(foldername(3)) ' ' num2str(foldername(4)) '.' num2str(foldername(5)) '_' num2str(APPForce) 'N'];
mkdir(foldername);
%eval(['cd ' foldername]);
path=pwd;
Folderpath = [path '\' foldername];

%%
delete 'Term.db';
delete 'Resorption.inp';
copyfile('Term-BAK.db','Term.db');
copyfile('Cell_Differentiation.m',Folderpath);
copyfile('preParameter.inp',Folderpath);
copyfile('2.Solve.inp',Folderpath);


OutHistoryType=[];
OutE=[];
StrainHist=[];
FluidHist=[];
E_10=[]; %1=ElemNum 2~11=E 12=AverageResult
P_10=[]; %1=ElemNum 2~11=E 12=AverageResult
Perm_10=[]; %1=ElemNum 2~11=E 12=AverageResult
resorpElem=[]; %�s��n�Q�l����Element

disp('-----Creating Folder Done -----')
disp('-----ParameterSetting Done -----')

for Iteration=1:teststep
    
    disp(['----------------Iteration : ' num2str(Iteration) ' Start----------------']);
    %% Pre-Parameter
    fid=fopen('Z_MATLAB_parameter.inp','w');
    fprintf(fid,'teststep = %d\n',Iteration);
    fprintf(fid,'APPForce = %d\n',APPForce);
    fclose(fid);
    %% ANSYS Part
    if(Iteration==1)
        %File 1.Presetting
        OperatingResult=system(['SET KMP_STACKSIZE=2048k & "' ansysPath '" -b -j "MATLAB_Differentiation" -i 1.Presetting.inp -o Output MATLAB_PART']);
        if OperatingResult~=8
            disp('File 1.Presetting .....Error');
            break;
        end
        disp('File 1.Presetting .....Done');
    end
    
    %File 2.Solve
    OperatingResult=system(['SET KMP_STACKSIZE=2048k & "' ansysPath '" -b -j "MATLAB_Differentiation" -i 2.Solve.inp -o Output MATLAB_PART']);
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
    FlowZ = RawData(:,end);
    Flow = RawData(:,7);
    EPX = RawData(:,8);
    EPY = RawData(:,9);
    EPZ = RawData(:,10);
    EPXY = RawData(:,11);
    EPYZ = RawData(:,12);
    EPXZ = RawData(:,13);
    Granu=Concentration<0.001;
    %Differentiation�����p��
    Octahedral_strain = sqrt((Prin1-Prin2).^2+(Prin2-Prin3).^2+(Prin3-Prin1).^2)./2; %Modeling Mechanical Signals on the Surface ofmCT and CAD Based Rapid Prototype ScaffoldModels to Predict (Early Stage) TissueDevelopment
%     Octahedral_strain = sqrt((EPX-EPY).^2+(EPY-EPZ).^2+(EPZ-EPX).^2+6.*((EPXY./2).^2+(EPYZ./2).^2+(EPXZ./2).^2)).*(2/3); %http://www.colorado.edu/MCEN/MCEN5023/chap_04.pdf
%     Octahedral_strain = sqrt((Prin1-Prin2).^2+(Prin2-Prin3).^2+(Prin3-Prin1).^2).*(2/3); %https://books.google.com.tw/books?id=PDF7INRk1tsC&pg=PA272&lpg=PA272&dq=Octahedral+shear+strain&source=bl&ots=UVKAMyIJEO&sig=mpQ1zDP9Mfav4lFSWFceG0YGzN0&hl=zh-TW&sa=X&ved=0ahUKEwiqzs22xtnPAhUMuo8KHQv6D3wQ6AEIfjAO#v=onepage&q=Octahedral%20shear%20strain&f=false
%     Octahedral_strain = abs(4.*Prin1.*Prin2.*Prin3).^(1/3); %Distortion stress, distortion strain and their physical concept
%     Octahedral_strain = sqrt(Prin1.^2-(Prin1.*Prin2)+(Prin2).^2).*sqrt(2/3); % A prediction of cell differentiation and proliferation within a collagen�Vglycosaminoglycan scaffold subjected to mechanical strain and perfusive fluid flow
%     Octahedral_strain = sqrt((Prin1-Prin2).^2+(Prin1-Prin3).^2+(Prin2-Prin3).^2)./sqrt(2); %A numerical model of the fracture healing process that describes tissue development and revascularisation

    
    Fluid_velocity = sqrt(abs((sqrt(FlowX.^2-FlowY.^2).^2 - FlowZ.^2)));
    %     Fluid_velocity = abs(Flow);
    S=(Octahedral_strain./a + Fluid_velocity./b);
    % Cell_type: 0=granulation tissue 1=mature bone 0.266>=S 2=immature bone 1>=S>0.266 3=cartilage 3>=S>1 4=fibrous connective tissue S>3
    DiffType = zeros(size(S));
    DiffType(S>3) = 4;
    DiffType(S>1 & S<=3) = 3;
    DiffType(S>0.266 & S<=1) = 2;
    DiffType(S>0.0103 & S<=0.266) = 1; %1=mature bone 0.266>=S
    DiffType(S<=0.0103) = 0; %0=resorption
    DiffType(Granu)=5;
    
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
        if DiffType(in)==5
            Concentration(in)=1;
        end
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
    
    %Smooth method = Get the average Young's modulus in 10 result (�Q�B����) 
    E_10(:,end)=Eaverage; %�x�}�Gsum(a)�O�C�D�M�Asum(a,2)�O��D�M ���B�D�X�Q�B�����é�b�̫�@���
    P_10(:,end)=Paverage; %�x�}�Gsum(a)�O�C�D�M�Asum(a,2)�O��D�M ���B�D�X�Q�B�����é�b�̫�@���
    Perm_10(:,end)=Permaverage; %�x�}�Gsum(a)�O�C�D�M�Asum(a,2)�O��D�M ���B�D�X�Q�B�����é�b�̫�@���
    
    
    UniMP=E_10(:,end-1+E_10_Option); %�s�����Young's modulus
    UniP=P_10(:,end-1+E_10_Option); %�s�����Poisson ratio
    UniPerm=Perm_10(:,end-1+E_10_Option); %�s�����Fluid permeability
    
    
    ProcessResult={}; %�Ĥ@��O���ưѼơA�᭱�d��O�ϥθӧ��ưѼƪ������s��
    
    
    for i=1:length(UniMP)
        ProcessResult{i}=[UniMP(i) UniP(i) UniPerm(i) Elenum(i)];
    end
    
    
    
    %% ��X���ưѼƤΤ����s��
    part=1;
    fid=fopen(['MP_Iteration_' num2str(part) '.inp'],'w');
    fprintf(fid,'resume,Term,db\n');
    fprintf(fid,'/prep7\n');
    fprintf(fid,'allsel\n');
    fprintf(fid,'/INPUT,Z_PreSettingClear,''inp'',,, 0\n');
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
    
    % ��X�Q�l����Element
    fid=fopen('Resorption.inp','w');
    rl=sum(DiffType==0);
    rlist=Elenum(DiffType==0);
    resorpElem=unique([resorpElem;rlist]);
    for r=1:length(resorpElem)
        fprintf(fid,'esel,u,,,%d\n',resorpElem(r));
    end
    fclose(fid);
       
    disp(['Material Properties Output ' num2str(part) '.....Done']);
    %% Plot Result
    fid=fopen([Folderpath '\Iteration_Result' '_' num2str(Iteration,'%04i') '.inp'],'w');
    fprintf(fid,'/prep7\n');
    fprintf(fid,'allsel\n');
    fprintf(fid,'et,6,185\n');
    fprintf(fid,'et,7,185\n');
    fprintf(fid,'et,8,185\n');
    fprintf(fid,'et,9,185\n');
    fprintf(fid,'et,10,185\n');
    fprintf(fid,'et,11,185\n');
    fprintf(fid,'/COLOR, 6, 0\n'); %Resorption
    fprintf(fid,'/COLOR, 7, 7\n'); %Mature Bone
    fprintf(fid,'/COLOR, 8, 9\n'); %Immature Bone
    fprintf(fid,'/COLOR, 9, 4\n'); %cartilage
    fprintf(fid,'/COLOR, 10, 2\n'); %fibtous
    fprintf(fid,'/COLOR, 11, 12\n'); %Granu
    for plotType=0:5
        if plotType~=0
            plotInd=(DiffType==plotType);
            plotL=sum(plotInd);
            plotElem=Elenum(plotInd);
        else
            plotElem=rlist;
            plotL=length(rlist);
            disp(['Resorp Num = ' num2str(plotL)]);
        end
        for in=2:plotL
        fprintf(fid,'emodif,%d,type,%d\n',plotElem(in),6+plotType);
        end
    end
    fclose(fid);
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
        mpr=system(['SET KMP_STACKSIZE=2048k & "' ansysPath '" -b -j "MATLAB_Differentiation" -i MP_Iteration_' num2str(MPpart) '.inp -o Output MP_input']);
        if mpr~=8
            error('MP input failed');
        end
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

% disp('Out put avi....')
% png2movie( Folderpath,'Iteration_Result__' );
end


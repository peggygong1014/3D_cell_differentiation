global ansysPath
%% °Ñ¼Æ³]¸m
addpath([pwd '\function']);
% preParameter;
% disp('Parameter definition is finished');
% 
% %% Build model
% delete('Term-BAK.db');
% OperatingResult=system(['SET KMP_STACKSIZE=2048k & "' ansysPath '" -b -j "Model" -dir "' ...
%     pwd  '\Model\" -i "' pwd '\Model\Z_Model_Build.inp" -o Output']);
% if OperatingResult~=8
%     disp('Model building is failed');
%     error('Model building is failed');
% end
% disp('Model have been built');
% %% Diffusion
% copyfile([pwd  '\Model\Term-BAK.db'],[pwd  '\Diffusion\'])
% copyfile([pwd  '\Model\Term-BAK.db'],pwd)

OperatingResult=system(['SET KMP_STACKSIZE=2048k & "' ansysPath '" -b -j "Term-BAK" -dir "' ...
    pwd  '\Diffusion\" -i "' pwd '\Diffusion\Diffusion.inp" -o Output']);

if OperatingResult~=8
    disp('Diffusion calculation is failed');
    error('Diffusion calculation is failed');
end

%% Algorithm
Cell_Differentiation;
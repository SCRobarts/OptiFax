%% cuda_compile

% clear all;
% cd('D:\Seb Robarts\OneDrive - Heriot-Watt University\Heriot Watt\Summer 2021\MATLABThings')
cd('C:\Users\Seb Robarts\Heriot-Watt University Team Dropbox\RES_EPS_McCracken_Lab\Seb\OptiFax\toolbox\lib\MEXlib\src')

mexcuda -dynamic -g  -O...
	-output OPOmexbatch512TPB ... % -output OPOmexbatchDynamic ...
	OPO_CUDA_BATCH_OP.cu ... % OPO_CUDA_DYNAMIC_BATCH_OP.cu ... 
	OPO_mex_BATCH.cpp...
	'-LC:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.2\lib' -lcufft

% copyfile OPOmexbatch.mexw64 ../OPOmexbatch.mexw64
copyfile OPOmexbatch512TPB.mexw64 ../OPOmexbatch512TPB.mexw64
% copyfile OPOmexbatchDynamic.mexw64 ../OPOmexbatchDynamic.mexw64


% 	-output OPO_mex_test ...
% 	OPO_TEST_CUDA_ADAPT.cu ...
% 	OPO_mex_test_adapt.cpp

% 	-DCUDA_API_PER_THREAD_DEFAULT_STREAM=1 ...

%  mexcuda '-LC:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.2\lib\x64' -lcufft -dynamic -O -g -DCUDA_API_PER_THREAD_DEFAULT_STREAM=1 -output OPO_mex_WIP OPO_CUDA_WIP_1_0.cu OPO_mex_WIP_1_0.cpp
%This file is downloaded from https://github.com/qunmingwang/Code-for-S2-fusion
%Copyright belongs to Dr. Qunming Wang

%%% from 3 dim to 2 dim.
function hyp_data2=D3_D2(hyp_data);
hyp_data2=shiftdim(hyp_data,2);
hyp_data2=hyp_data2(:,:);
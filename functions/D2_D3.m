%This file is downloaded from https://github.com/qunmingwang/Code-for-S2-fusion
%Copyright belongs to Dr. Qunming Wang

%%% from 2 dim to 3 dim.
function hyp_data3=D2_D3(hyp_data2,line,column);
[dim,nn]=size(hyp_data2);
hyp_data3=zeros(line,column,dim);
for i=1:dim
    hyp_data3(:,:,i)=reshape(hyp_data2(i,:),line,column);
end
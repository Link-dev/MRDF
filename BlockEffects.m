function data = BlockEffects(prediction, origin, w, numS)
%%% mitigate the block effects by neighborhood similar pixels
%
%The similar pixels are belong to the same class with the center pixel
%The weight of similar pixel is determined by spatial euclidean distance
%
%This code uses parallel computing
%
%parameters:
%    prediciton:        the result with block effects
%    origin:            the origin fine image with few bands
%    w:                 the half window size, if 25, the window size is 25*2+1=51
%    numS:              number of similar pixels
%
%
% sample:
%    data = BlockEffects(prediction, orgin, w, numS)

if nargin<3 || isempty(w)
    w=25;
end
if nargin<4 || isempty(numS)
    numS=20;         
end

[~, ~, nbp] = size(prediction);
[ns, nl, nb] = size(origin);
data = zeros(ns,nl, nbp);

parfor i = 1:ns
    for j = 1:nl

        ai = max(1, i-w); bi = min(ns, i+w);
        aj = max(1,j-w);  bj = min(nl, j+w);
        blockSize = (bi-ai+1)*(bj-aj+1);
        ci = i-ai+1; cj = j-aj+1;
        
        fineTmp = origin(ai:bi, aj:bj, :);
        num_similar_pixel1 = min(numS,blockSize);
        tmp2 = sqrt(mean((reshape(fineTmp, blockSize, nb) - reshape(fineTmp(ci,cj,:), 1, nb)).^2,2));
        
        [~,I] = sort(tmp2);
        tmp2(I(num_similar_pixel1+1:blockSize)) = Inf; 
        tmp2 = reshape(tmp2, (bi-ai+1), (bj-aj+1));
        [col,row] = find(tmp2 ~= Inf);
        
        Distance = 1 + sqrt((col-ci).^2 + (row-cj).^2) / w;
        Weighted = 1./Distance/sum(1./Distance);
        
        tmp = reshape(prediction(ai:bi, aj:bj, :), blockSize, nbp);
        tmp = tmp(I(1:numS), :);
        
        data(i,j,:) = tmp'*Weighted;
    end
end

end
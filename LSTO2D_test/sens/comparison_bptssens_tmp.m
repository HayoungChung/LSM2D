
% format long
LOG = zeros(100,2);
for cnt = 0:99

A = textread(sprintf('cant_bpts_%d.txt',cnt)); 
B = textread(sprintf('lsq_bpts_%d.txt',cnt)); 
% % B = B(:,[3,4,1,2]);
% Bdot = zeros(size(B));
% for ii = 1:4
% tmp = reshape(B(:,ii),80,[])';
% Bdot(:,ii) = tmp(:);
% end

LOG(cnt+1,:) = [cnt, max(abs(A(:)./B(:)))];


end
% format short
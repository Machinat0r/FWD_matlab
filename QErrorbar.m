function QErrorbar(varargin)
% 
%------written by Wending Fu, Jun.26.2022 in Beijing------------
t = varargin{1};
Q = varargin{2};
resY = varargin{3};

%% plot
plot(t,Q,'k','LineWidth',1.5);hold on;

%% error
for i = 1:length(Q)
posY(i,:) = resY{i};negY(i,:) = resY{i};
posY(i,posY(i,:)<=Q(i)) = nan;
negY(i,negY(i,:)>Q(i)) = nan;
% posY(i,:) = resY{i}(resY{i}>Q(i));
% negY(i,:) = resY{i}(resY{i}<=Q(i));
end
maxY = max(posY,[],2);
minY = min(negY,[],2);

%% vertical line
c_eval("line([t(?),t(?)],[minY(?),maxY(?)],'color','k');hold on;",1:length(Q))
%% scatter max & min
errorbar(t,Q,minY-Q,maxY-Q,'k');
% scatter(t,maxY,'kv','MarkerFaceColor','w');hold on;
% scatter(t,minY,'k^','MarkerFaceColor','w');hold on;
%% scatter other
posY(posY == maxY) = nan;
negY(negY == minY) = nan;
for i = 1:length(Q)
    scatter(repmat(t(i),1,size(posY(i,:),2)),posY(i,:),'k+');hold on;
    scatter(repmat(t(i),1,size(posY(i,:),2)),negY(i,:),'k+');hold on;
end

end
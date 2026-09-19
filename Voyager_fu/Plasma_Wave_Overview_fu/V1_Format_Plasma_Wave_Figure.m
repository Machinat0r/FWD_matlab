function V1_Format_Plasma_Wave_Figure(f)
% 仅整理标题与正常图例位置；不改变任何科学数据或坐标范围。
titleText='Voyager 1 | 1990-01-01 to 2025-06-30 | P1 0.57-1.78 MeV';
items=findall(f,'-property','String');
for k=1:numel(items)
    if ~isgraphics(items(k)), continue, end
    value=items(k).String;
    if (ischar(value) || (isstring(value) && isscalar(value))) && strcmp(string(value),titleText)
        delete(items(k));
    end
end
annotation(f,'textbox',[.085 .953 .83 .036],'String',titleText, ...
    'EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',15,'FontWeight','bold','Tag','OverviewTitle');
lg=findall(f,'Type','legend');
assert(numel(lg)==1,'Expected the single EPO/QTN density legend.');
lg.Units='normalized';lg.Position=[.128 .373 .14 .029];
drawnow;
end

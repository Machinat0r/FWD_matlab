function MMS_orbit_layout_20260924(fig,h)
% 将MMS_orbit原生的8个坐标轴横向排版；不改变数据、模型或坐标定义。
set(fig,'Renderer','painters','PaperPositionMode','auto','Position',[0 0 1600 1050]);
for ip=1:8
    col=mod(ip-1,4);row=floor((ip-1)/4);
    set(h(ip),'Units','normalized','Position',[0.055+0.245*col 0.56-0.47*row 0.185 0.36]);
end
set(h([1 3 4]),'XLim',[-35 20]);
set(h(1:7),'FontSize',18);
set(h(8),'Position',[0.78 0.09 0.21 0.36]);
tx=findall(h(8),'Type','text');
for it=1:numel(tx)
    str=string(tx(it).String);
    if any(str==["MMS1","MMS2","MMS3","MMS4"])
        sc=str2double(extractAfter(str,'MMS'));
        tx(it).Position=[0.04+0.55*mod(sc-1,2) 1.0-0.13*floor((sc-1)/2) 0];
        tx(it).FontSize=17;
    elseif contains(str,'MMS configuration')
        tx(it).Position=[0 0.71 0];tx(it).FontSize=17;
    elseif contains(str,'IMF')
        tx(it).Position=[0 0.43 0];tx(it).FontSize=17;
        tx(it).String=strrep(strrep(char(str),',By=',',\newline By='),',Bz=',',\newline Bz=');
    end
end
ln=findall(h(8),'Type','line');
for il=1:numel(ln)
    x=ln(il).XData;
    if numel(x)==1
        sc=round(x/0.27)+1;
        if sc>=1&&sc<=4
            ln(il).XData=0.55*mod(sc-1,2);
            ln(il).YData=1.0-0.13*floor((sc-1)/2);
        end
    end
end
end

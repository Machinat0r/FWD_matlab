clear;clc
ParentDir='C:\Matlab\bin\新建文件夹\fwd\Sundry\zwq\GOLD\';
DirStruct = dir(ParentDir);

 for z=3:length(DirStruct)
     filename = [ParentDir DirStruct(z).name];
%      filename = [ParentDir, 'GOLD_L2_TDISK_2020_', repmat('0',1,3-ceil(log10(z+0.1))),num2str(z), '_v03_r01_c01.mat'];
%      filename = strrep(filename,' ','');

 
     load(filename,'tdisk');
% load('C:\matgold\GOLD_L2_TDISK_2020_230_v03_r01_c01.mat');
%     h = fspecial('average',3);
b=size(tdisk,3);
 for i=1:b
    p=tdisk(:,:,i);
%     figure(1); pcolor(p);
%     colorbar;
%  caxis([300,600]);
    ppp=p(2:45,2:51);
    ppe=p(2:45,3:52);
    pps=p(3:46,2:51);
    ppw=p(2:45,1:50);
    ppn=p(1:44,2:51);
    ppen=p(1:44,3:52);
    ppwn=p(1:44,1:50);
    ppes=p(3:46,3:52);
    ppws=p(3:46,1:50);
    for m=1:44
        for n=1:50
        if ((~isnan(ppp(m,n))))
            juzheng=[ppe(m,n) pps(m,n) ppw(m,n) ppn(m,n),ppp(m,n),ppen(m,n),ppwn(m,n),ppes(m,n),ppws(m,n)];
            ppp(m,n)=nanmean(juzheng);
        end
        end
    end
    p(2:45,2:51)=ppp;
%     figure(2); pcolor(p);colorbar;
%     caxis([300,600]);
    tdisk(:,:,i)=p;
 end
 Tdisk = tdisk;
 save(filename, 'Tdisk', '-append');
 end
    
%     hahah
%     ppp=p(2:45,2:51);
%     ppe=p(2:45,3:52);
%     pps=p(3:46,2:51);
%     ppw=p(2:45,1:50);
%     ppn=p(1:44,2:51);
%     ppen=p(1:44,3:52);
%     ppwn=p(1:44,1:50);
%     ppes=p(3:46,3:52);
%     ppws=p(3:46,1:50);
%     for m=1:44
%         for n=1:50
%         if ((~isnan(ppp(m,n))) )
%             juzheng=[ppe(m,n) pps(m,n) ppw(m,n) ppn(m,n),ppp(m,n),ppen(m,n),ppwn(m,n),ppes(m,n),ppws(m,n)];
%             ppp(m,n)=nanmean(juzheng);
%         end
%         end
%     end
%     p(2:45,2:51)=ppp;
%     figure(3); pcolor(p);
%     colorbar;
%  caxis([300,600]);
%  haha
%     
    
    
%     pp=p(2:45,2:51);
%     pe=p(2:45,3:52);
%     ps=p(3:46,2:51);
%     pw=p(2:45,1:50);
%     pn=p(1:44,2:51);
%     for m=1:44
%         for n=1:50
%         if (isnan(pp(m,n))) && ((~isnan(pe(m,n)))||(~isnan(ps(m,n)))||(~isnan(pw(m,n)))||(~isnan(pn(m,n))))
%             juzheng=[pe(m,n) ps(m,n) pw(m,n) pn(m,n)];
%             pp(m,n)=nanmean(juzheng);
%         end
%     end
%     end
%     bbb=pp(1,:);
%     pp=[bbb;pp];
%     aaa=pp(:,1);
%     pp=[aaa pp];
%     bbb=pp(1,:);
%     pp=[pp;bbb];
%     aaa=pp(:,1);
%     pp=[pp aaa];
%     pp_=pp;
%     
%     pp = imfilter(pp,h,'replicate');
%     
%     ppp=pp(2:45,2:51);
%     ppe=pp(2:45,3:52);
%     pps=pp(3:46,2:51);
%     ppw=pp(2:45,1:50);
%     ppn=pp(1:44,2:51);
%         for m=1:44
%         for n=1:50
%         if (isnan(ppp(m,n))) && ((~isnan(ppe(m,n)))||(~isnan(pps(m,n)))||(~isnan(ppw(m,n)))||(~isnan(ppn(m,n))))
%             juzheng=[ppe(m,n) pps(m,n) ppw(m,n) ppn(m,n)];
%             ppp(m,n)=nanmean(juzheng);
%         end
%     end
%         end
%      bbb=ppp(1,:);
%     ppp=[bbb;ppp];
%     aaa=ppp(:,1);
%     ppp=[aaa ppp];
%     bbb=ppp(1,:);
%     ppp=[ppp;bbb];
%     aaa=ppp(:,1);
%     ppp=[ppp aaa];
%     %将矩阵东南西北各平移一下，生成四个矩阵
%     %h = fspecial('log',[3 3],0.5) ;
% 
%     p_smooth = imfilter(ppp,h,'replicate');
% % end
% figure(1);pcolor(p)
% 
% figure(2);pcolor(p_smooth)

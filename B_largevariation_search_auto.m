clear;clc
Date = '2017-07-01/2017-07-05'; %不会爬取endDate的数据，如若要2017.7.4则将终止日期设为2017.7.5
splitDate = regexp(Date,'/','split');
ic = 1:4;

filenames = SDCFilenames(Date,ic,'inst','fgm','drm','brst');
% filenames2 = SDCFilenames(Date,ic,'inst','fpi','drm','brst','dpt','des-moms,dis-moms,des-dist,dis-dist');
% filenames3 = SDCFilenames(Date,ic,'inst','scm','drm','brst','dpt','scb');
% filenames_srvy = SDCFilenames(Date,ic,'inst','fgm','drm','srvy'); %为了知道坐标
% filenames = [filenames1,filenames2,filenames3];

global OutputDir
OutputDir = ['C:\Matlab\MMS\',splitDate{1},'To',splitDate{2},'\'];
DataPath = [OutputDir,'Data\'];
mkdir(DataPath)
SDCFilesDownload(filenames,DataPath)
SDCDataMove(DataPath,'C:\Matlab\MMS\')

expr = '_[0-9]+\_v';
NameTags = regexp(filenames,expr,'match');
date = {};
for i = 1:length(NameTags)
    date{end+1} = NameTags{i}{1}(2:9);
end
date = unique(date);

t_start=clock;
for tempd = 1:length(date)
%----------------------------------------------------%
%读取四颗卫星磁场数据列表
ic=1:4;
c_eval(['fileFolder?=fullfile(''C:\Matlab\MMS\mms?\fgm\brst\l2\',date{tempd}(1:4),...
    '\',date{tempd}(5:6),'\',date{tempd}(7:8),'''',');'],ic);
c_eval('dirOutput?=dir(fullfile(fileFolder?,''*''));',ic);
c_eval('name?=transpose({dirOutput?.name});',ic);


try %防止某颗卫星某天完全没数据而在寻找共有磁场时报错
    c_eval('name?(1:2,:)=[];',ic);
    c_eval('name?=char(name?);',ic);
    c_eval('filesize?=size(name?);',ic);
catch
    writematrix([date{tempd},'缺失一颗或多颗卫星数据'],[OutputDir,'errorlog.txt'],'WriteMode','append','Encoding','UTF-8')
    continue 
end



%----------------------------------------------------%
%找寻四颗卫星共有的磁场数据
c_eval('name?_temp=name?(:,18:31);',ic);
c_eval('name?_temp=str2num(name?_temp);',ic);
intersection_MMS1_4=intersect(intersect(intersect(name1_temp,name2_temp),name3_temp),name4_temp);
[~,~,id1]=intersect(intersection_MMS1_4,name1_temp);
[~,~,id2]=intersect(intersection_MMS1_4,name2_temp);
[~,~,id3]=intersect(intersection_MMS1_4,name3_temp);
[~,~,id4]=intersect(intersection_MMS1_4,name4_temp);
c_eval('nameB?=name?(id?,:);',ic);


% c_eval('fileFolderi?=fullfile(''F:\DF数据\fpidis\mms?'');',ic);
% c_eval('dirOutputi?=dir(fullfile(fileFolderi?,''*''));',ic);
% c_eval('namei?=transpose({dirOutputi?.name});',ic);
% c_eval('namei?(1:2,:)=[];',ic);
% c_eval('namei?=char(namei?);',ic);
% c_eval('filesizei?=size(namei?);',ic);
% 
% scnume=3;%有没有第四颗卫星的电子数据，有scnmme=4,否则是3；
% if scnume==3
%     ic=1:3;
% else
%     ic=1:4;
% end
% 
% c_eval('fileFoldere?=fullfile(''F:\DF数据\fpides\mms?'');',ic);
% c_eval('dirOutpute?=dir(fullfile(fileFoldere?,''*''));',ic);
% c_eval('namee?=transpose({dirOutpute?.name});',ic);
% c_eval('namee?(1:2,:)=[];',ic);
% c_eval('namee?=char(namee?);',ic);
% c_eval('filesizee?=size(namee?);',ic);


ic=1:4;

%----------------------------------------------------%
% 读取数据并判断
num=0;%满足条件的个数
eventlist=cell(60000,4);
filesize=length(id1);
for i=1:filesize
    % 读取数据
    ic=1:4;
    c_eval(['route?=''C:\Matlab\MMS\mms?\fgm\brst\l2\',date{tempd}(1:4),...
    '\',date{tempd}(5:6),'\',date{tempd}(7:8),'\''',';'],ic);
    clear filename1 filename2 filename3 filename4;
    c_eval('filename?=strcat(route?,nameB?(i,:));',ic);
    clear fpiFile1 fpiFile2 fpiFile3 fpiFile4 B1 B2 B3 B4; 
    c_eval('fpiFile?=dataobj(filename?);',ic);
    c_eval('B?=getmat(fpiFile?,''mms?_fgm_b_gsm_brst_l2'');',ic);
    
    c_eval('marktime?=regexp(nameB?(i,:),''_'',''start'');',ic);
%     c_eval('namei?1=strcat(''mms?_fpi_brst_l2_dis-moms_'',nameB?(i,marktime?(4)+1:marktime?(5)-1));',ic);
%     c_eval('namee?1=strcat(''mms?_fpi_brst_l2_des-moms_'',nameB?(i,marktime?(4)+1:marktime?(5)-1));',ic);
%     
%     for j=1:filesizei1(1,1)
%            if ~strcmp(namei11,namei1(j,1:length(namei11)))==0
%                namei12=namei1(j,:);
%            end
%     end
%     
%     for j=1:filesizei2(1,1)
%            if ~strcmp(namei21,namei2(j,1:length(namei21)))==0
%                namei22=namei2(j,:);
%            end
%     end
%     
%     for j=1:filesizei3(1,1)
%            if ~strcmp(namei31,namei3(j,1:length(namei31)))==0
%                namei32=namei3(j,:);
%            end
%     end
%     
%     for j=1:filesizei4(1,1)
%            if ~strcmp(namei41,namei4(j,1:length(namei41)))==0
%                namei42=namei4(j,:);
%            end
%     end
%     
%     c_eval('routei?=''F:\DF数据\fpidis\mms?\'';',ic);
%     clear filenamei1 filenamei2 filenamei3 filenamei4;
%     c_eval('filenamei?=strcat(routei?,namei?2(1,:));',ic);
%     
%     clear fpiFilei1 fpiFilei2 fpiFilei3 fpiFilei4 Vi1 Vi2 Vi3 Vi4; 
%     c_eval('fpiFilei?=dataobj(filenamei?);',ic);
%     c_eval('Vigse?=getmat(fpiFilei?,''mms?_dis_bulkv_gse_brst'');',ic);
%     c_eval('Vi?=irf_gse2gsm(Vigse?);',ic);
%     
%     for j=1:filesizee1(1,1)
%            if ~strcmp(namee11,namee1(j,1:length(namee11)))==0
%                namee12=namee1(j,:);
%            end
%     end
%     
%     for j=1:filesizee2(1,1)
%            if ~strcmp(namee21,namee2(j,1:length(namee21)))==0
%                namee22=namee2(j,:);
%            end
%     end
%     
%     for j=1:filesizee3(1,1)
%            if ~strcmp(namee31,namee3(j,1:length(namee31)))==0
%                namee32=namee3(j,:);
%            end
%     end
%     
%     if scnume==4
%         for j=1:filesizee4(1,1)
%            if ~strcmp(namee41,namee4(j,1:length(namee41)))==0
%                namee42=namee4(j,:);
%            end
%         end
%     end
%     
%     if scnume==3
%         ic=1:3;
%     else
%         ic=1:4;
%     end
%     
%     c_eval('routee?=''F:\DF数据\fpides\mms?\'';',ic);
%     clear filenamee1 filenamee2 filenamee3 filenamee4;
%     c_eval('filenamee?=strcat(routee?,namee?2(1,:));',ic);
%     
%     clear fpiFilee1 fpiFilee2 fpiFilee3 fpiFilee4 Ve1 Ve2 Ve3 Ve4; 
%     c_eval('fpiFilee?=dataobj(filenamee?);',ic);
%     c_eval('Vegse?=getmat(fpiFilee?,''mms?_des_bulkv_gse_brst'');',ic);
%     c_eval('Ve?=irf_gse2gsm(Vegse?);',ic);
%     
%     if scnume==3
%         Ve4=[Ve3(:,1) 0*Ve3(:,2:end)];
%     end
%     
    
       
    ic=2:4;
    c_eval('B?=irf_resamp(B?,B1);',ic);
    %四颗卫星磁场大小差值
    deltB12=[B1(:,1) abs(B1(:,5)-B2(:,5))];
    deltB13=[B1(:,1) abs(B1(:,5)-B3(:,5))];
    deltB14=[B1(:,1) abs(B1(:,5)-B4(:,5))];
    deltB23=[B1(:,1) abs(B2(:,5)-B3(:,5))];
    deltB24=[B1(:,1) abs(B2(:,5)-B4(:,5))];
    deltB34=[B1(:,1) abs(B3(:,5)-B4(:,5))];
    %四颗卫星磁场夹角
    angle12=[B1(:,1) acosd(irf_dot(B1(:,1:4),B2(:,1:4),1)./(B1(:,5).*B2(:,5)))];
    angle13=[B1(:,1) acosd(irf_dot(B1(:,1:4),B3(:,1:4),1)./(B1(:,5).*B3(:,5)))];
    angle14=[B1(:,1) acosd(irf_dot(B1(:,1:4),B4(:,1:4),1)./(B1(:,5).*B4(:,5)))];
    angle23=[B1(:,1) acosd(irf_dot(B2(:,1:4),B3(:,1:4),1)./(B2(:,5).*B3(:,5)))];
    angle24=[B1(:,1) acosd(irf_dot(B2(:,1:4),B4(:,1:4),1)./(B2(:,5).*B4(:,5)))];
    angle34=[B1(:,1) acosd(irf_dot(B3(:,1:4),B4(:,1:4),1)./(B3(:,5).*B4(:,5)))];
    
    
    widln=1920;   %128 Hz，窗口长度15s
    indsta=1; indend=indsta+widln;
    dataln=length(B1(:,1));
    for indsta=1:widln:dataln
        
        clear id id_thr id_diff id_con flag;
        %在这个窗口里，差值大于0.5nT或者角度差10°的点
        deltB_threshold=1;
        angleB_threshold=10;
        if ((dataln-indend)>=widln)
            
            indend=indsta+widln;
            
%             deltB_threshold=1;
%             angleB_threshold=10;
            id=find((deltB12(indsta:indend,2)>=deltB_threshold)|(deltB13(indsta:indend,2)>=deltB_threshold)|...
                (deltB14(indsta:indend,2)>=deltB_threshold)|(deltB23(indsta:indend,2)>=deltB_threshold)|...
                (deltB24(indsta:indend,2)>=deltB_threshold)|(deltB34(indsta:indend,2)>=deltB_threshold)|...
                (angle12(indsta:indend,2)>=angleB_threshold)|(angle13(indsta:indend,2)>=angleB_threshold)|...
                (angle14(indsta:indend,2)>=angleB_threshold)|(angle23(indsta:indend,2)>=angleB_threshold)|...
                (angle24(indsta:indend,2)>=angleB_threshold)|(angle34(indsta:indend,2)>=angleB_threshold));
            id_thr=id+(indsta-1)*zeros(length(id),1);
            id_diff=diff(id_thr);
            id_con=find(id_diff==1);%如果出现连续的大于阈值的点，差值为1
            flag=sum(id_diff(id_con));%1的个数
            if flag>=64
                timesta=B1(indsta,1);
                timeend=B1(indend,1); 
                ic=1:4;
                c_eval('Btint?=irf_tlim(B?,[timesta timeend]);',ic);
                plot_B(Btint1,Btint2,Btint3,Btint4,timesta,timeend);
                num=num+1;
                timestautc=irf_time(timesta,'epoch>utc');
                timeendutc=irf_time(timeend,'epoch>utc');
                datastautc=irf_time(B1(1,1),'epoch>utc');
                dataendutc=irf_time(B1(end,1),'epoch>utc');
                eventlist{num,1}=num;eventlist{num,2}=timestautc;eventlist{num,3}=timeendutc;
                eventlist{num,4}=datastautc;eventlist{num,5}=dataendutc;
            end
             
        else
            widln1=dataln-indend;
            indend=dataln;
             %在这个窗口里，差值大于0.5nT或者角度差10°的点
%             deltB_threshold=0.5;
%             angleB_threshold=10;
            id=find((deltB12(indsta:indend,2)>=deltB_threshold)|(deltB13(indsta:indend,2)>=deltB_threshold)|...
                (deltB14(indsta:indend,2)>=deltB_threshold)|(deltB23(indsta:indend,2)>=deltB_threshold)|...
                (deltB24(indsta:indend,2)>=deltB_threshold)|(deltB34(indsta:indend,2)>=deltB_threshold)|...
                (angle12(indsta:indend,2)>=angleB_threshold)|(angle13(indsta:indend,2)>=angleB_threshold)|...
                (angle14(indsta:indend,2)>=angleB_threshold)|(angle23(indsta:indend,2)>=angleB_threshold)|...
                (angle24(indsta:indend,2)>=angleB_threshold)|(angle34(indsta:indend,2)>=angleB_threshold));
            id_thr=id+(indsta-1)*zeros(length(id),1);
            id_diff=diff(id_thr);
            id_con=find(id_diff==1);%如果出现连续的大于阈值的点，差值为1
            flag=sum(id_diff(id_con));%1的个数
            if flag>=widln1/20
                timesta=B1(indsta,1);
                timeend=B1(indend,1); 
                ic=1:4;
                c_eval('Btint?=irf_tlim(B?,[timesta timeend]);',ic);
                plot_B(Btint1,Btint2,Btint3,Btint4,timesta,timeend);
                num=num+1;
                timestautc=irf_time(timesta,'epoch>utc');
                timeendutc=irf_time(timeend,'epoch>utc');
                datastautc=irf_time(B1(1,1),'epoch>utc');
                dataendutc=irf_time(B1(end,1),'epoch>utc');
                eventlist{num,1}=num;eventlist{num,2}=timestautc;eventlist{num,3}=timeendutc;
                eventlist{num,4}=datastautc;eventlist{num,5}=dataendutc;
            end
                
        end
        
    end
       
end
end
% eventlist(1,1)=1111;eventlist(1,2)=2222;eventlist(1,3)=3333;
% eventlist(1,4)=4444;eventlist(1,5)=5555;
try
    xlswrite([OutputDir,'evenlist'],eventlist);
    t_end=clock;
    etime(t_end,t_start)
catch
    '';
end

function plot_B(dataB1,dataB2,dataB3,dataB4,datatimesta,datatimeend)
Tsta=datatimesta;
Tend=datatimeend;
tint=[Tsta Tend];
Tstautc=irf_time(Tsta,'epoch>utc');
Tendutc=irf_time(Tend,'epoch>utc');
%% Init figure 1
n_subplots=4;
i_subplot=1;
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLineLineWidth', 0.5);
fn=figure(1);clf;
set(gcf,'PaperUnits','centimeters');
xSize = 18; ySize = 30; coef=floor(min(800/xSize,800/ySize));
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[10 10 xSize*coef ySize*coef]);

%% B plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([dataB1(:,1) dataB1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([dataB2(:,1) dataB2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([dataB3(:,1) dataB3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([dataB4(:,1) dataB4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([dataB1(:,1) dataB1(:,2)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('Bx [nT]','fontsize',10);
set(gca,'ColorOrder',[[0 0 0];[1 0 0];[0 1 0];[0 0 1]]);
irf_legend(gca,{'MMS1','MMS2','MMS3','MMS4'},[0.1 0.12]);
% irf_legend(gca,'a',[0.99 0.98],'color','k','fontsize',10);
%% B plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([dataB1(:,1) dataB1(:,3)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([dataB2(:,1) dataB2(:,3)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([dataB3(:,1) dataB3(:,3)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([dataB4(:,1) dataB4(:,3)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([dataB1(:,1) dataB1(:,3)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('By [nT]','fontsize',10);
% irf_legend(gca,'b',[0.99 0.98],'color','k','fontsize',10);
%% B plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([dataB1(:,1) dataB1(:,4)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([dataB2(:,1) dataB2(:,4)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([dataB3(:,1) dataB3(:,4)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([dataB4(:,1) dataB4(:,4)], 'color','b', 'Linewidth',0.75); hold on;
irf_plot([dataB1(:,1) dataB1(:,4)*0],'k--', 'Linewidth',0.75);hold off;
grid off;
ylabel('Bz [nT]','fontsize',10);
% irf_legend(gca,'c',[0.99 0.98],'color','k','fontsize',10);
%% B plot
h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
irf_plot([dataB1(:,1) dataB1(:,5)], 'color','k', 'Linewidth',0.75); hold on;
irf_plot([dataB2(:,1) dataB2(:,5)], 'color','r', 'Linewidth',0.75); hold on;
irf_plot([dataB3(:,1) dataB3(:,5)], 'color','g', 'Linewidth',0.75); hold on;
irf_plot([dataB4(:,1) dataB4(:,5)], 'color','b', 'Linewidth',0.75); hold on;
% irf_plot([dataB1(:,1) dataB1(:,5)*0],'k--', 'Linewidth',0.75);
hold off;
grid off;
ylabel('|B| [nT]','fontsize',10);
% irf_legend(gca,'d',[0.99 0.98],'color','k','fontsize',10);
% %% Vi plot
% h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% irf_plot([dataVi1(:,1) dataVi1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([dataVi2(:,1) dataVi2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([dataVi3(:,1) dataVi3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([dataVi4(:,1) dataVi4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
% irf_plot([dataVi1(:,1) dataVi1(:,2)*0],'k--', 'Linewidth',0.75);hold off;
% grid off;
% ylabel('Vix [km/s]','fontsize',10);
% %% Vi plot
% h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% irf_plot([dataVi1(:,1) dataVi1(:,3)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([dataVi2(:,1) dataVi2(:,3)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([dataVi3(:,1) dataVi3(:,3)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([dataVi4(:,1) dataVi4(:,3)], 'color','b', 'Linewidth',0.75); hold on;
% irf_plot([dataVi1(:,1) dataVi1(:,3)*0],'k--', 'Linewidth',0.75);hold off;
% grid off;
% ylabel('Viy [km/s]','fontsize',10);
% %% Vi plot
% h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% irf_plot([dataVi1(:,1) dataVi1(:,4)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([dataVi2(:,1) dataVi2(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([dataVi3(:,1) dataVi3(:,4)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([dataVi4(:,1) dataVi4(:,4)], 'color','b', 'Linewidth',0.75); hold on;
% irf_plot([dataVi1(:,1) dataVi1(:,4)*0],'k--', 'Linewidth',0.75);hold off;
% grid off;
% ylabel('Viz [km/s]','fontsize',10);
% %% Ve plot
% h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% irf_plot([dataVe1(:,1) dataVe1(:,2)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([dataVe2(:,1) dataVe2(:,2)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([dataVe3(:,1) dataVe3(:,2)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([dataVe4(:,1) dataVe4(:,2)], 'color','b', 'Linewidth',0.75); hold on;
% irf_plot([dataVe1(:,1) dataVe1(:,2)*0],'k--', 'Linewidth',0.75);hold off;
% grid off;
% ylabel('Vex [km/s]','fontsize',10);
% %% Ve plot
% h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% irf_plot([dataVe1(:,1) dataVe1(:,3)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([dataVe2(:,1) dataVe2(:,3)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([dataVe3(:,1) dataVe3(:,3)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([dataVe4(:,1) dataVe4(:,3)], 'color','b', 'Linewidth',0.75); hold on;
% irf_plot([dataVe1(:,1) dataVe1(:,3)*0],'k--', 'Linewidth',0.75);hold off;
% grid off;
% ylabel('Vey [km/s]','fontsize',10);
% %% Ve plot
% h(i_subplot)=irf_subplot(n_subplots,1,-i_subplot);i_subplot=i_subplot+1;
% irf_plot([dataVe1(:,1) dataVe1(:,4)], 'color','k', 'Linewidth',0.75); hold on;
% irf_plot([dataVe2(:,1) dataVe2(:,4)], 'color','r', 'Linewidth',0.75); hold on;
% irf_plot([dataVe3(:,1) dataVe3(:,4)], 'color','g', 'Linewidth',0.75); hold on;
% irf_plot([dataVe4(:,1) dataVe4(:,4)], 'color','b', 'Linewidth',0.75); hold on;
% irf_plot([dataVe1(:,1) dataVe1(:,4)*0],'k--', 'Linewidth',0.75);hold off;
% grid off;
% ylabel('Vez [km/s]','fontsize',10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
legendtick='abcdefghijkl';
for kk=1:i_subplot-1
    irf_legend(h(kk),legendtick(kk),[0.99 0.98],'color','k','fontsize',10);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(h(1:end),'fontsize',10);
irf_adjust_panel_position;
irf_plot_ylabels_align(h(1:end));
irf_zoom(tint,'x',h(1:end));
%% 
set(gcf,'render','painters');
% pause(1)
global OutputDir
FigPath = [OutputDir,'Fig\'];
mkdir(FigPath);
% set(gcf,'paperpositionmode','auto')
figname=[ FigPath,'Events_B','_' ,Tstautc(1:4), Tstautc(6:7), Tstautc(9:10), '-', Tstautc(12:13), ...
    Tstautc(15:16), Tstautc(18:19), '_' ,Tendutc(12:13), Tendutc(15:16), Tendutc(18:19)];
print(gcf, '-dpng', [figname '.png']);
close;
end

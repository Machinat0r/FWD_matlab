


function c=c_rapid_subspin_pad_plot_bf2004(varargin)
% c_rapid_subspin_pad_plot_bf2004 plot the RAPID differential flux under subspin resolution
%
%Input parameter (3 in total, downloaded from CAA website and then constructed into the object format)
%1. C?_CP_RAP_EPADEX?
%2. C?_CP_RAP_EPITCH
%3. E_channel=1;             channel 1 is 40.7 keV; channel 2 is 68.1keV
%4. tint;             time interval you want to plot
%
%see also c_rapid_pad_plot_bf2004
%---------------------------------------
%Example: 
%c_rapid_subspin_pad_plot_bf2004(C1_CP_RAP_EPADEX2, C1_CP_RAP_EPITCH, 1, tint);
%---------------------------------------

%----created by Huishan Fu at IRFU (2012-06-30)----


[ax,args,nargs] = axescheck(varargin{:});
RAPflux_obj=args{1};
RAPpitang_obj=args{2};
E_channel=args{3};
tint=args{4};


    varsFlux=RAPflux_obj.Variables;
    Flux_str=varsFlux{3,1};  Azm_str=varsFlux{6,1};  Time_str=varsFlux{2,1}; 
    varsPitang=RAPpitang_obj.Variables;
    Pitang_str=varsPitang{3,1};

    RapFlux=getmat(RAPflux_obj, Flux_str);
    RapPitang=getmat(RAPpitang_obj,Pitang_str);
    RapAzmang=getmat(RAPflux_obj, Azm_str);
    Raptime=getmat(RAPflux_obj, Time_str);
    
    aaa=RapFlux.data;
    bbb=RapPitang.data;
    ttt=Raptime(:,1);
    
    Polang=[10 30 50 70 90 110 130 150 170];
    
    Flux_subspin=zeros(length(ttt)*16,9);
    Pitang_subspin=zeros(length(ttt)*16,9);
    Time_subspin=linspace(ttt(1),ttt(end),length(ttt)*16);
    
    
    for kk=1:length(ttt)
        Flux_each=squeeze(aaa(kk,E_channel,:,:));
        Pitang_each=squeeze(bbb(kk,:,:));
        
        ind_range=[(16*(kk-1)+1):(16*kk)];
        Flux_subspin(ind_range,:)=Flux_each;
        Pitang_subspin(ind_range,:)=Pitang_each;
    end
    
    
    indsta=find(Time_subspin>=tint(1)); indsta=indsta(1);
    indend=find(Time_subspin>=tint(2)); indend=indend(1);
    
    Time_subspin=Time_subspin(1,indsta:indend);
    Pitang_subspin=Pitang_subspin(indsta:indend,:);
    Flux_subspin=Flux_subspin(indsta:indend,:);
    
    Xaxis=[indsta:indend];
    
    specFlux=pcolor(gca, Xaxis, linspace(0,180,10), ([Flux_subspin Flux_subspin(:,end)])');
    set(specFlux,'EdgeColor','none');
    hold on;
    colorbar;
    colormap(jet);

 
    [C PAcontour]=contour(gca, Xaxis,Polang,Pitang_subspin',[30 90 150]); hold off;
    clabel(C, PAcontour, [30 90 150], 'FontSize', 8, 'Rotation',0);
    set(PAcontour, 'color','r', 'Linewidth',0.7);

    
end




    
   
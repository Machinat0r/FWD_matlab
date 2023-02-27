


function c=c_rapid_subspin_flux_bf2004(varargin)
% c_rapid_subspin_flux_bf2004 plot the RAPID differential flux under subspin resolution
%
%Input parameter (2 in total, downloaded from CAA website and then constructed into the object format)
%1. C?_CP_RAP_EPADEX?
%2. E_channel=1;             channel 1 is 40.7 keV; channel 2 is 68.1keV
%
%see also c_rapid_pad_plot_bf2004
%---------------------------------------
%Example: 
%c_rapid_subspin_flux_bf2004(C1_CP_RAP_EPADEX2, 1);
%---------------------------------------

%----created by Huishan Fu at IRFU (2012-06-30)----


[ax,args,nargs] = axescheck(varargin{:});
RAPflux_obj=args{1};
E_channel=args{2};


    varsFlux=RAPflux_obj.Variables;
    Flux_str=varsFlux{3,1};  Azm_str=varsFlux{6,1};  Time_str=varsFlux{2,1}; 

    RapFlux=getmat(RAPflux_obj, Flux_str);
    RapAzmang=getmat(RAPflux_obj, Azm_str);
    Raptime=getmat(RAPflux_obj, Time_str);
    
    aaa=RapFlux.data;
    ttt=Raptime(:,1);
    
    Polang=[10 30 50 70 90 110 130 150 170];
    
    Flux_subspin=zeros(length(ttt)*16,9);
    Time_subspin=linspace(ttt(1),ttt(end),length(ttt)*16);
    
    
    for kk=1:length(ttt)
        Flux_each=squeeze(aaa(kk,E_channel,:,:));
        
        ind_range=[(16*(kk-1)+1):(16*kk)];
        Flux_subspin(ind_range,:)=Flux_each;
    end
    
    
    Flux_subspin_interg=nanmean(Flux_subspin,2);
    irf_plot(gca, [Time_subspin' Flux_subspin_interg]);


end




    
   
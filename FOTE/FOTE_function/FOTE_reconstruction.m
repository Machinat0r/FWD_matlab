function [h,FOTE_res] = FOTE_reconstruction(dB_null, Null_loc, focus_points, BoxWid, view_angles, plottype, cmin, cmax)
fig = figure();
h=[];
h(1)=axes('position',[0.1 0.1 0.8 0.8]);
Siz = size(focus_points);
for is = 1:Siz(1)
    Xgrid = focus_points(is,1); Ygrid = focus_points(is,2); Zgrid = focus_points(is,3);
    %inverse trace
    ii = 1;
    Xprev=Xgrid;  Yprev=Ygrid;  Zprev=Zgrid;
    while abs(Xprev)<=BoxWid && abs(Yprev)<=BoxWid && abs(Zprev)<=BoxWid
        Xcurt=Xprev; 
        Ycurt=Yprev; 
        Zcurt=Zprev;
        Br = dB_null * ([Xcurt Ycurt Zcurt])';
        
        Bxcurt=Br(1);
        Bycurt=Br(2);
        Bzcurt=Br(3);
        Bmcurt=sqrt(Br(1).^2+Br(2).^2+Br(3).^2);
        
%         step=Bmcurt;
        step=1;
        stepvec=[Bxcurt Bycurt Bzcurt]/norm([Bxcurt Bycurt Bzcurt])*step;
        Xprev=Xcurt-stepvec(1);
        Yprev=Ycurt-stepvec(2);
        Zprev=Zcurt-stepvec(3);
        
        Xline(ii)=Xcurt;
        Yline(ii)=Ycurt;
        Zline(ii)=Zcurt;
        Bmline(ii)=Bmcurt;
        ii=ii+1;
    end
    
    Xline_inv(is,1:length(Xline)) = Xline;
    Yline_inv(is,1:length(Xline)) = Yline;
    Zline_inv(is,1:length(Xline)) = Zline;
    Bmline_inv(is,1:length(Xline)) = Bmline;
    
    clear Xline
    clear Yline
    clear Zline
    clear Bmline
    
    %positive trace
    ii=1;
    Xnext=Xgrid;  Ynext=Ygrid;  Znext=Zgrid;
    
    while abs(Xnext)<=BoxWid && abs(Ynext)<=BoxWid && abs(Znext)<=BoxWid
        Xcurt=Xnext;
        Ycurt=Ynext;
        Zcurt=Znext;
        
        Br = dB_null * ([Xcurt Ycurt Zcurt])';
        Bxcurt=Br(1);
        Bycurt=Br(2);
        Bzcurt=Br(3);
        Bmcurt=sqrt(Br(1).^2+Br(2).^2+Br(3).^2);
        
%         step=Bmcurt;
        step=1;
        stepvec=[Bxcurt Bycurt Bzcurt]/norm([Bxcurt Bycurt Bzcurt])*step;
        Xnext=Xcurt+stepvec(1);
        Ynext=Ycurt+stepvec(2);
        Znext=Zcurt+stepvec(3);
        
        Xline(ii)=Xcurt;
        Yline(ii)=Ycurt;
        Zline(ii)=Zcurt;
        Bmline(ii)=Bmcurt;
        ii=ii+1;
    end
    
    Xline_pos(is,1:length(Xline)) = Xline;
    Yline_pos(is,1:length(Xline)) = Yline;
    Zline_pos(is,1:length(Xline)) = Zline;
    Bmline_pos(is,1:length(Xline)) = Bmline;
    
    clear Xline
    clear Yline
    clear Zline
    clear Bmline
end
pa=1;na=1;arrP=0.5;ap=15;
for nk = 1:Siz(1)
    ns0=find(Xline_inv(nk,:)~=0); a=size(ns0);a=a(2);
    if a~=0
        ns0 = ns0(length(ns0));
%         if Zline_inv(nk,ns0)<0
            switch plottype
                case('line')
                    plot3(gca, Xline_inv(nk,1:ns0), Yline_inv(nk,1:ns0), Zline_inv(nk,1:ns0),'b'); hold on;
                otherwise
                    cline(Xline_inv(nk,1:ns0), Yline_inv(nk,1:ns0), Zline_inv(nk,1:ns0), Bmline_inv(nk,1:ns0), cmin, cmax, jet);hold on;
            end
            % arrow3matrix_preva(pa,1:3)=[Xline_preva(nk,fix(arrP*ns0)),Yline_preva(nk,fix(arrP*ns0)),Zline_preva(nk,fix(arrP*ns0))];
            % arrow3matrix_preva(pa,4:6)=[Xline_preva(nk,fix(arrP*ns0)-ap),Yline_preva(nk,fix(arrP*ns0)-ap),Zline_preva(nk,fix(arrP*ns0)-ap)];
            pa = pa+1;
%         end
        ns0=find(Xline_pos(nk,:)~=0); ns0 = ns0(length(ns0));
%         if Zline_pos(nk,ns0)<0
            switch plottype
                case('line')
                    plot3(gca, Xline_pos(nk,1:ns0), Yline_pos(nk,1:ns0), Zline_pos(nk,1:ns0),'r'); hold on;
                otherwise
                    cline(Xline_pos(nk,1:ns0), Yline_pos(nk,1:ns0), Zline_pos(nk,1:ns0), Bmline_pos(nk,1:ns0), cmin, cmax, jet);hold on;
            end
            % arrow3matrix_nexta(na,1:3)=[Xline_nexta(nk,fix(arrP*ns0)),Yline_nexta(nk,fix(arrP*ns0)),Zline_nexta(nk,fix(arrP*ns0))];
            % arrow3matrix_nexta(na,4:6)=[Xline_nexta(nk,fix(arrP*ns0)+ap),Yline_nexta(nk,fix(arrP*ns0)+ap),Zline_nexta(nk,fix(arrP*ns0)+ap)];
            na = na+1;
%         end
    end
end
%% SC PLOT
Rsc1 = Null_loc.R_sc1; Rsc2 = Null_loc.R_sc2;
Rsc3 = Null_loc.R_sc3; Rsc4 = Null_loc.R_sc4;
plot3(gca, [Rsc1(1) Rsc2(1) Rsc3(1) Rsc4(1)], [Rsc1(2) Rsc2(2) Rsc3(2) Rsc4(2)], ...
[Rsc1(3) Rsc2(3) Rsc3(3) Rsc4(3)], 'k', 'Linewidth',1); hold on;
plot3(gca, [Rsc2(1) Rsc4(1) Rsc1(1) Rsc3(1)], [Rsc2(2) Rsc4(2) Rsc1(2) Rsc3(2)], ...
[Rsc2(3) Rsc4(3) Rsc1(3) Rsc3(3)], 'k', 'Linewidth',1); hold on;
plot3(gca, [Rsc1(1)], [Rsc1(2)],[Rsc1(3)], 'ks', 'Linewidth',1, ...
'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',12); hold on;
plot3(gca, [Rsc2(1)], [Rsc2(2)],[Rsc2(3)], 'rs', 'Linewidth',1, ...
'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',12); hold on;
plot3(gca, [Rsc3(1)], [Rsc3(2)],[Rsc3(3)], 'gs', 'Linewidth',1, ...
'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',12); hold on;
plot3(gca, [Rsc4(1)], [Rsc4(2)],[Rsc4(3)], 'bs', 'Linewidth',1, ...
'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',12); hold on;

set(gca,'view',view_angles)
box on
Xline_inv(Xline_inv == 0) = nan;
Yline_inv(Yline_inv == 0) = nan;
Zline_inv(Zline_inv == 0) = nan;
Xline_pos(Xline_pos == 0) = nan;
Yline_pos(Yline_pos == 0) = nan;
Zline_pos(Zline_pos == 0) = nan;
inv_p = struct('X',Xline_inv,'Y',Yline_inv,'Z',Zline_inv); pos_p = struct('X',Xline_pos,'Y',Yline_pos,'Z',Zline_pos);
FOTE_res.inv_p = inv_p; FOTE_res.pos_p = pos_p;
end
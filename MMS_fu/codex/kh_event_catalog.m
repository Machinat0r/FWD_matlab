function events = kh_event_catalog()
%KH_EVENT_CATALOG Published MMS magnetopause Kelvin-Helmholtz intervals.
%
% The base catalog is the 45 clear KHI encounters in Rice et al. (2022),
% DOI: 10.1029/2021JA029685.  Rows are promoted or broadened when later
% case studies provide stronger rolled-up-vortex evidence.  Additional
% non-overlapping intervals come from published case studies and from the
% MMS rows of the Kavosi et al. (2023) supplementary event list.
%
% Confidence values intentionally distinguish mature rolled-up vortices
% from clear KHI/KHW encounters and probable/linear cases.

tz = 'UTC';
parse = @(s) datetime(s,'InputFormat','yyyy-MM-dd''T''HH:mm:ss','TimeZone',tz);

% id, start, end, side, source, doi
rice = {
 'R01','2015-09-08T09:00:00','2015-09-08T11:50:00','dusk','Rice et al. 2022; Wilder et al. 2023','10.1029/2021JA029685;10.1029/2023JA031583'
 'R02','2015-09-15T10:45:00','2015-09-15T14:45:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R03','2015-10-11T10:30:00','2015-10-11T11:00:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R04','2015-10-15T06:00:00','2015-10-15T07:00:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R05','2015-10-17T16:00:00','2015-10-17T16:28:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R06','2015-10-18T15:00:00','2015-10-18T15:25:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R07','2015-12-22T22:15:00','2015-12-22T22:50:00','dawn','Rice et al. 2022','10.1029/2021JA029685'
 'R08','2016-01-11T20:52:00','2016-01-11T21:10:00','dawn','Rice et al. 2022','10.1029/2021JA029685'
 'R09','2016-01-19T19:57:00','2016-01-19T20:35:00','dawn','Rice et al. 2022','10.1029/2021JA029685'
 'R10','2016-02-05T18:55:00','2016-02-05T19:30:00','dawn','Rice et al. 2022','10.1029/2021JA029685'
 'R11','2016-02-07T03:45:00','2016-02-07T04:40:00','dawn','Rice et al. 2022','10.1029/2021JA029685'
 'R12','2016-02-18T19:30:00','2016-02-18T20:40:00','dawn','Rice et al. 2022','10.1029/2021JA029685'
 'R13','2016-02-25T18:55:00','2016-02-25T20:05:00','dawn','Rice et al. 2022; Nykyri et al. 2021','10.1029/2021JA029685;10.1029/2019JA027698'
 'R14','2016-09-26T14:15:00','2016-09-26T15:25:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R15','2016-09-27T18:00:00','2016-09-27T20:10:00','dusk','Rice et al. 2022; Tang et al. 2018; Hou et al. 2026; Hwang et al. 2022; Wilder et al. 2023','10.1029/2021JA029685;10.5194/angeo-36-879-2018;10.3389/fspas.2026.1837800;10.3389/fspas.2022.895514;10.1029/2023JA031583'
 'R16','2016-10-04T18:20:00','2016-10-04T19:30:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R17','2016-10-10T14:40:00','2016-10-10T15:40:00','dusk','Rice et al. 2022; Guo et al. 2026','10.1029/2021JA029685;10.26464/epp2026042'
 'R18','2016-10-24T10:50:00','2016-10-24T11:20:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R19','2016-11-04T11:45:00','2016-11-04T13:00:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R20','2017-05-03T02:00:00','2017-05-03T04:30:00','dawn','Rice et al. 2022','10.1029/2021JA029685'
 'R21','2017-05-08T13:00:00','2017-05-08T14:50:00','dawn','Rice et al. 2022','10.1029/2021JA029685'
 'R22','2017-05-11T12:00:00','2017-05-11T14:30:00','dawn','Rice et al. 2022','10.1029/2021JA029685'
 'R23','2017-05-11T14:44:03','2017-05-11T16:14:03','dawn','Rice et al. 2022; Wilder et al. 2023','10.1029/2021JA029685;10.1029/2023JA031583'
 'R24','2017-05-19T23:58:00','2017-05-20T01:45:00','dawn','Rice et al. 2022','10.1029/2021JA029685'
 'R25','2017-05-20T02:00:00','2017-05-20T04:30:00','dawn','Rice et al. 2022','10.1029/2021JA029685'
 'R26','2017-09-20T22:32:23','2017-09-20T22:57:03','dusk','Rice et al. 2022; Wilder et al. 2023','10.1029/2021JA029685;10.1029/2023JA031583'
 'R27','2017-09-26T16:35:00','2017-09-26T19:07:00','dusk','Rice et al. 2022; Wilder et al. 2023','10.1029/2021JA029685;10.1029/2023JA031583'
 'R28','2017-10-16T14:30:00','2017-10-16T15:20:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R29','2017-10-30T19:05:00','2017-10-30T19:40:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R30','2017-11-02T17:25:00','2017-11-02T18:15:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R31','2018-05-03T00:15:00','2018-05-03T00:50:00','dawn','Rice et al. 2022','10.1029/2021JA029685'
 'R32','2018-09-18T15:50:00','2018-09-18T16:15:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R33','2018-09-24T14:10:00','2018-09-24T17:30:00','dusk','Rice et al. 2022; Wilder et al. 2023','10.1029/2021JA029685;10.1029/2023JA031583'
 'R34','2018-10-02T23:45:00','2018-10-03T00:20:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R35','2018-10-04T17:10:00','2018-10-04T18:00:00','dusk','Rice et al. 2022; Kavosi et al. 2023','10.1029/2021JA029685;10.1038/s41467-023-37485-x'
 'R36','2019-04-13T07:45:00','2019-04-13T08:15:00','dawn','Rice et al. 2022','10.1029/2021JA029685'
 'R37','2019-06-03T23:05:00','2019-06-04T00:20:00','dawn','Rice et al. 2022','10.1029/2021JA029685'
 'R38','2019-09-25T13:45:00','2019-09-26T02:30:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R39','2019-10-02T08:15:00','2019-10-02T11:00:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R40','2019-10-02T16:00:00','2019-10-02T17:20:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R41','2019-10-02T21:40:00','2019-10-02T22:05:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R42','2019-10-06T14:00:00','2019-10-06T18:30:00','dusk','Rice et al. 2022; Wilder et al. 2023','10.1029/2021JA029685;10.1029/2023JA031583'
 'R43','2019-10-15T19:00:00','2019-10-15T20:15:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R44','2019-10-22T22:00:00','2019-10-22T22:20:00','dusk','Rice et al. 2022','10.1029/2021JA029685'
 'R45','2019-11-12T20:28:00','2019-11-12T20:45:00','dusk','Rice et al. 2022; Guo et al. 2026','10.1029/2021JA029685;10.26464/epp2026042'
};

n = size(rice,1);
events = table(string(rice(:,1)),parse(string(rice(:,2))),parse(string(rice(:,3))), ...
    strings(n,1),strings(n,1),string(rice(:,4)),string(rice(:,5)),string(rice(:,6)), ...
    repmat("B_clear_KHI",n,1),repmat("clear KHI/KHW; rolled-up state not established for every row",n,1), ...
    'VariableNames',{'CatalogID','StartUTC','EndUTC','PreferredStartText','PreferredEndText','Side','Source','DOI','Tier','Confidence'});

% Promote rows supported by rolled-up-vortex case studies.
promote = ["R15","R23","R26","R27","R33","R42","R45"];
events.Tier(ismember(events.CatalogID,promote)) = "A_rolled_up";
events.Confidence(ismember(events.CatalogID,promote)) = "published case study reports nonlinear/rolled-up KH structures";
events.Confidence(events.CatalogID=="R01") = "clear KHI; mature nonlinear roll-up not established (Wilder 2023)";
events.Tier(events.CatalogID=="R13") = "C_probable";
events.Confidence(events.CatalogID=="R13") = "KHI-consistent high-latitude boundary layer; alternative pressure forcing discussed";
events.Tier(events.CatalogID=="R17") = "C_linear";
events.Confidence(events.CatalogID=="R17") = "KH structures reported as linear or not fully rolled up";

events.PreferredStartText(events.CatalogID=="R01") = "2015-09-08T10:27:00";
events.PreferredEndText(events.CatalogID=="R01")   = "2015-09-08T10:47:00";
events.PreferredStartText(events.CatalogID=="R15") = "2016-09-27T19:50:00";
events.PreferredEndText(events.CatalogID=="R15")   = "2016-09-27T20:06:00";
events.PreferredStartText(events.CatalogID=="R23") = "2017-05-11T15:44:00";
events.PreferredEndText(events.CatalogID=="R23")   = "2017-05-11T16:04:00";
events.PreferredStartText(events.CatalogID=="R26") = "2017-09-20T22:35:00";
events.PreferredEndText(events.CatalogID=="R26")   = "2017-09-20T22:55:00";
events.PreferredStartText(events.CatalogID=="R27") = "2017-09-26T17:50:00";
events.PreferredEndText(events.CatalogID=="R27")   = "2017-09-26T18:10:00";
events.PreferredStartText(events.CatalogID=="R33") = "2018-09-24T15:50:00";
events.PreferredEndText(events.CatalogID=="R33")   = "2018-09-24T16:10:00";
events.PreferredStartText(events.CatalogID=="R42") = "2019-10-06T16:05:00";
events.PreferredEndText(events.CatalogID=="R42")   = "2019-10-06T16:25:00";

% Extra non-overlapping, peer-reviewed MMS intervals.
extra = {
 'X01','2015-09-10T19:40:00','2015-09-10T21:10:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X02','2015-09-14T15:50:00','2015-09-14T16:10:00','','','dusk','Hwang et al. 2022','10.3389/fspas.2022.895514','B_KHV_table','coordinated-observation table labels magnetopause KHVs; engineering window around listed time'
 'X03','2015-10-01T17:59:00','2015-10-01T18:09:00','','','dusk','Hwang et al. 2022','10.3389/fspas.2022.895514','A_rolled_up','case-study magnetopause KHVs'
 'X04','2015-10-08T15:55:00','2015-10-08T16:15:00','','','dusk','Hwang et al. 2022','10.3389/fspas.2022.895514','B_KHV_table','coordinated-observation table labels magnetopause KHVs; engineering window around listed time'
 'X05','2016-02-06T18:33:00','2016-02-06T20:15:00','2016-02-06T19:47:00','2016-02-06T19:54:00','dawn','Hwang et al. 2022','10.3389/fspas.2022.895514','A_rolled_up','multiple-spacecraft confirmation of dawnside KH vortices'
 'X06','2016-02-07T18:42:00','2016-02-07T19:02:00','','','dawn','Hwang et al. 2022','10.3389/fspas.2022.895514','B_KHV_table','coordinated-observation table labels magnetopause KHVs; engineering window around listed time'
 'X07','2016-02-18T17:39:00','2016-02-18T17:59:00','','','dawn','Hwang et al. 2022','10.3389/fspas.2022.895514','B_KHV_table','coordinated-observation table labels magnetopause KHVs; engineering window around listed time'
 'X08','2016-10-23T16:10:00','2016-10-23T17:00:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X09','2016-11-01T18:35:00','2016-11-01T19:10:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X10','2016-11-20T17:00:00','2016-11-20T17:25:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X11','2017-02-17T15:28:00','2017-02-17T15:48:00','','','unknown','Hwang et al. 2022','10.3389/fspas.2022.895514','B_KHV_table','coordinated-observation table labels magnetopause KHVs; engineering window around listed time'
 'X12','2017-05-05T19:20:00','2017-05-05T23:20:00','2017-05-05T20:05:00','2017-05-05T20:25:00','dawn','Kieokaew et al. 2020; Hwang et al. 2020','10.1029/2019JA027527;10.1029/2019JA027665','A_rolled_up','nonlinear rolled-up dawnside KH vortices'
 'X13','2017-05-14T05:55:00','2017-05-14T06:10:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X14','2017-05-28T22:30:00','2017-05-28T22:55:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X15','2017-05-29T00:30:00','2017-05-29T01:50:00','2017-05-29T01:00:00','2017-05-29T01:20:00','dawn','Wilder et al. 2023','10.1029/2023JA031583','A_rolled_up','low-density high-speed signature confirms rolled-up KHI'
 'X16','2017-05-29T04:30:00','2017-05-29T07:00:00','2017-05-29T05:20:00','2017-05-29T05:40:00','dawn','Lu et al. 2019; Kavosi et al. 2023','10.3847/1538-4357/ab0e76;10.1038/s41467-023-37485-x','C_linear','KHW/LDHS interval; not reported as highly nonlinear'
 'X17','2017-09-03T15:35:00','2017-09-03T16:00:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X18','2017-09-06T00:00:00','2017-09-06T00:20:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X19','2017-09-06T01:20:00','2017-09-06T01:50:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X20','2017-09-16T18:20:00','2017-09-16T18:45:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X21','2017-09-16T20:00:00','2017-09-16T20:45:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X22','2017-09-19T14:15:00','2017-09-19T16:00:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X23','2017-09-23T15:30:00','2017-09-23T16:30:00','2017-09-23T15:33:00','2017-09-23T15:53:00','dusk','Blasl et al. 2022','10.1063/5.0067370','A_rolled_up','southward-IMF KH waves and a rolled-up vortex'
 'X24','2018-04-18T20:30:00','2018-04-18T21:00:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X25','2018-04-18T21:20:00','2018-04-18T22:00:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X26','2018-04-30T03:55:00','2018-04-30T04:20:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X27','2018-04-30T04:40:00','2018-04-30T06:20:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X28','2018-05-03T01:40:00','2018-05-03T02:10:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X29','2018-05-03T04:10:00','2018-05-03T04:40:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X30','2018-05-08T21:00:00','2018-05-08T22:20:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X31','2018-05-16T00:30:00','2018-05-16T01:45:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X32','2018-05-20T22:20:00','2018-05-20T23:50:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X33','2018-05-28T20:40:00','2018-05-28T21:45:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X34','2018-09-30T00:40:00','2018-09-30T01:40:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X35','2018-10-08T22:25:00','2018-10-08T23:20:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X36','2018-10-09T01:30:00','2018-10-09T03:00:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X37','2018-10-14T17:40:00','2018-10-14T18:30:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X38','2018-10-14T19:00:00','2018-10-14T21:00:00','','','unknown','Kavosi et al. 2023','10.1038/s41467-023-37485-x','B_clear_KHW','catalogued MMS KHW'
 'X39','2019-05-20T21:41:23','2019-05-20T22:13:23','2019-05-20T21:50:00','2019-05-20T22:10:00','dawn','Wilder et al. 2023','10.1029/2023JA031583','A_rolled_up','low-density high-speed signature confirms rolled-up KHI'
 'X40','2019-10-12T14:30:03','2019-10-12T15:33:03','2019-10-12T14:55:00','2019-10-12T15:15:00','dawn','Wilder et al. 2023','10.1029/2023JA031583','A_rolled_up','low-density high-speed signature confirms rolled-up KHI'
 'X41','2020-06-26T00:12:00','2020-06-26T01:22:00','2020-06-26T00:40:00','2020-06-26T01:00:00','dawn','Li et al. 2023; Settino et al. 2026','10.1029/2023GL105539;10.1029/2025JA034952','A_rolled_up','nonlinear KHW with repeated magnetopause crossings and reconnection'
 'X42','2021-05-23T04:00:00','2021-05-23T05:40:00','2021-05-23T04:40:00','2021-05-23T05:00:00','dusk','Wilder et al. 2023','10.1029/2023JA031583','A_rolled_up','low-density high-speed signature confirms rolled-up KHI'
 'X43','2021-06-15T07:00:00','2021-06-15T09:00:00','2021-06-15T07:50:00','2021-06-15T08:10:00','dawn','Wilder et al. 2023','10.1029/2023JA031583','A_rolled_up','low-density high-speed signature confirms rolled-up KHI'
 'X44','2022-04-14T11:32:00','2022-04-14T12:40:00','2022-04-14T12:20:00','2022-04-14T12:40:00','dawn','Gurram et al. 2026','10.3389/fspas.2026.1788081','A_rolled_up','storm-time early nonlinear/rolled-up KH waves'
};

m = size(extra,1);
extraT = table(string(extra(:,1)),parse(string(extra(:,2))),parse(string(extra(:,3))), ...
    string(extra(:,4)),string(extra(:,5)),string(extra(:,6)),string(extra(:,7)), ...
    string(extra(:,8)),string(extra(:,9)),string(extra(:,10)), ...
    'VariableNames',events.Properties.VariableNames);
events = [events; extraT];
events = sortrows(events,'StartUTC');

% Parse optional preferred windows after concatenation.
events.PreferredStartUTC = NaT(height(events),1,'TimeZone',tz);
events.PreferredEndUTC = NaT(height(events),1,'TimeZone',tz);
for k = 1:height(events)
    if strlength(events.PreferredStartText(k)) > 0
        events.PreferredStartUTC(k) = parse(events.PreferredStartText(k));
    end
    if strlength(events.PreferredEndText(k)) > 0
        events.PreferredEndUTC(k) = parse(events.PreferredEndText(k));
    end
end
events.PreferredStartText = [];
events.PreferredEndText = [];

events.EventID = "KH" + compose('%03d',(1:height(events))');
events = movevars(events,'EventID','Before','CatalogID');
events.DurationMinutes = minutes(events.EndUTC-events.StartUTC);
events = movevars(events,'DurationMinutes','After','EndUTC');
end

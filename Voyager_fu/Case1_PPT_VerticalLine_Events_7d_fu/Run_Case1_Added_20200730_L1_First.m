function result = Run_Case1_Added_20200730_L1_First(cadence)
%Run_Case1_Added_20200730_L1_First Replot only this supplementary V1 case.
%   No argument generates both overviews and the five-time hourly PAD.
%   Existing scientific settings are supplied by the shared runner.
%   Modified: 2026-09-04
if nargin < 1, cadence = 'both'; end
eventDate = datetime(2020, 7, 30, 'TimeZone', 'UTC');
selectionAudit = struct;
selectionAudit.Request = 'Additional event circled in the user-provided 2020 annual overview.';
selectionAudit.DateBasis = 'Original COHO P1 local maximum:2020-07-30 07:00 UTC,3.587; nearby magnetic magnitude reaches0.537 nT onJuly29 18:00 andJuly30 03:00 UTC.';
selectionAudit.InputImage = 'C:\Users\ADMINI~1\AppData\Local\Temp\codex-clipboard-8f3ef5cb-d952-48ab-a10d-fb779faba15b.png';
result = Run_Case1_Additional_L1_First(eventDate, cadence, selectionAudit);
end

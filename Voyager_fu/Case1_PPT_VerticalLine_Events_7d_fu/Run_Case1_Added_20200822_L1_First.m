function result = Run_Case1_Added_20200822_L1_First(cadence)
%Run_Case1_Added_20200822_L1_First Replot only this supplementary V1 case.
%   No argument generates both overviews and the five-time hourly PAD.
%   Existing scientific settings are supplied by the shared runner.
%   Modified: 2026-09-04
if nargin < 1, cadence = 'both'; end
eventDate = datetime(2020, 8, 22, 'TimeZone', 'UTC');
selectionAudit = struct;
selectionAudit.Request = 'Additional event at the first (leftmost) blue vertical line in the 2020 overview.';
selectionAudit.DateBasis = 'Raster line is nearAug22 and has no exact timestamp. Original COHO magnetic magnitude reaches0.530 nT onAug22, with local P1 flux maximum3.023 onAug22 17:00 UTC. This anchors the event day; the five-time PAD peak is independently selected inside its seven-day window.';
selectionAudit.InputImage = 'C:\Users\ADMINI~1\AppData\Local\Temp\codex-clipboard-5ce6f56e-d3fc-4b36-8cde-659986b4d7fb.png';
result = Run_Case1_Additional_L1_First(eventDate, cadence, selectionAudit);
end

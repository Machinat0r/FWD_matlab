function MMS_event_overview_20260923_controller(spacecraft)
if nargin<1,spacecraft=1:4;end
% 隐藏MATLAB进程：只对候选CDF已全部下载的事件/卫星作图。
codeRoot=fileparts(mfilename('fullpath'));addpath(codeRoot);
audit='Z:\SPART-WORK\Data\MMS\derived\event_overview_20260923';
M=jsondecode(fileread(fullfile(audit,'manifest.json')));
events=jsondecode(fileread(fullfile(codeRoot,'MMS_event_overview_20260923_events.json')));
attempts=zeros(15,4);
while true
    pending=0;progress=false;
    for ie=1:15
        for ic=spacecraft
            base=sprintf('%s_%s_MMS%d_overview',events(ie).id,strrep(events(ie).start(1:10),'-',''),ic);
            done=fullfile(audit,[base '.json']);
            if isfile(done),continue;end
            pending=pending+1;
            if attempts(ie,ic)>=3,continue;end
            cov=M.candidate_coverage(strcmp({M.candidate_coverage.event},events(ie).id) & ...
                startsWith({M.candidate_coverage.product},sprintf('mms%d_',ic)));
            wanted={};
            for k=1:numel(cov)
                x=cov(k).candidate_files;
                if ischar(x),x={x};end
                wanted=[wanted;x(:)]; %#ok<AGROW>
            end
            ready=true;
            for k=1:numel(wanted)
                j=find(strcmp({M.files.file_name},wanted{k}),1);
                if isempty(j),ready=false;break;end
                f=dir(M.files(j).path);
                if isempty(f)||f.bytes~=M.files(j).file_size,ready=false;break;end
            end
            if ~ready,continue;end
            fprintf('READY event %d MMS%d (%d source files)\n',ie,ic,numel(wanted));
            attempts(ie,ic)=attempts(ie,ic)+1;
            MMS_event_overview_20260923(ie,ic);
            progress=true;
        end
    end
    count=numel(dir(fullfile(audit,'EV*_overview.json')));
    fid=fopen(fullfile(audit,sprintf('plot_progress_sc%d.json',spacecraft(1))),'w','n','UTF-8');
    fprintf(fid,'%s',jsonencode(struct('completed',count,'expected',60,'pending',pending,'attempts',attempts)));
    fclose(fid);
    if pending==0,break;end
    if isfile(fullfile(audit,'download_complete.json')) && ~progress
        fprintf('Download finished, plotting still has %d pending jobs; inspect missing/read errors.\n',pending);
        break
    end
    pause(15);
end
fprintf('PLOT_CONTROLLER_FINISHED\n');
end

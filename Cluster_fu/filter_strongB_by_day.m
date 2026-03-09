function filter_strongB_by_day(inFile, outFile)
%FILTER_STRONGB_BY_DAY  按天筛选 Strong B 记录并输出到新文件
%
% 用法：
%   filter_strongB_by_day('input.txt', 'filtered.txt')
%
% 规则：
%   1) 同一天必须同时存在 maxB1~maxB4（至少各一条），否则删掉当天
%   2) 同一天 maxB1~4 中任一值 > 1000，删掉当天
%   3) 同一天 maxB1~4 的时间戳中，最大-最小 > 1小时，删掉当天
%   4) 若 maxB1~maxB4 均不>100（即都 <=100），删掉当天
%
% 说明：
%   - 若同一天某个 maxBk 有多条记录，默认取该天该 k 的“最大值”那条作为代表
%   - 规则2/3/4 均基于这 4 条代表记录判断

    if nargin < 2
        error('用法：filter_strongB_by_day(inFile, outFile)');
    end

    lines = readlines(inFile, "EmptyLineRule", "skip");

    % 正则：日期(yyyymmdd) 时间(HH:MM:SS(.mmm可选)) k(1-4) 数值(支持科学计数法)
    pattern = 'Strong B at:\s*(\d{8})\s*(\d{2}:\d{2}:\d{2}(?:\.\d{1,3})?)\s*maxB([1-4])\s*=\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)';

    dateDay = datetime.empty(0,1);
    dtFull  = datetime.empty(0,1);
    bIdx    = zeros(0,1);
    bVal    = zeros(0,1);
    rawLine = strings(0,1);

    % 解析每行
    for i = 1:numel(lines)
        s = strtrim(lines(i));
        tk = regexp(s, pattern, 'tokens', 'once');
        if isempty(tk)
            continue; % 不匹配的行忽略
        end

        dStr = tk{1};                    % yyyymmdd
        tStr = normalize_time_str(tk{2}); % HH:MM:SS.SSS (补齐到3位毫秒)
        k    = str2double(tk{3});        % 1..4
        v    = str2double(tk{4});

        day = datetime(dStr, "InputFormat", "yyyyMMdd");
        dt  = datetime(dStr + " " + tStr, "InputFormat", "yyyyMMdd HH:mm:ss.SSS");

        dateDay(end+1,1) = day; %#ok<AGROW>
        dtFull(end+1,1)  = dt;  %#ok<AGROW>
        bIdx(end+1,1)    = k;   %#ok<AGROW>
        bVal(end+1,1)    = v;   %#ok<AGROW>
        rawLine(end+1,1) = s;   %#ok<AGROW>
    end

    if isempty(rawLine)
        warning('未解析到任何匹配的记录。');
        writelines(strings(0,1), outFile);
        return;
    end

    % 按天分组
    [G, days] = findgroups(dateDay);

    outLines = strings(0,1);

    for g = 1:numel(days)
        idxDay = find(G == g);

        selIdx = nan(1,4); % 每天选出的 maxB1..4 的代表记录索引

        % 规则1：必须 1..4 都存在；并为每个 k 选择“值最大”的那条
        ok = true;
        for k = 1:4
            idxK = idxDay(bIdx(idxDay) == k);
            if isempty(idxK)
                ok = false;
                break;
            end
            [~, imax] = max(bVal(idxK));
            selIdx(k) = idxK(imax);
        end
        if ~ok
            continue; % 删掉当天
        end

        % 规则2：任一 > 1000 删掉当天
        vSel = bVal(selIdx);
        if any(vSel > 1000)
            continue;
        end

        % 规则4：若 maxB1~4 均不>100（即都 <= 100）删掉当天
        if all(vSel <= 100)
            continue;
        end

        % 规则3：时间跨度 > 1h 删掉当天（用选出的 4 条代表记录）
        tSel = dtFull(selIdx);
        if max(tSel) - min(tSel) > hours(1)
            continue;
        end

        % 通过筛选：输出这一天选出的4条记录（按 B1..B4 顺序）
        outLines = [outLines; rawLine(selIdx(:))]; %#ok<AGROW>
    end

    writelines(outLines, outFile);
end

function tOut = normalize_time_str(tIn)
% 将 "HH:MM:SS" 或 "HH:MM:SS.m" / ".mm" / ".mmm" 统一为 "HH:MM:SS.SSS"
    tIn = string(tIn);
    if ~contains(tIn, ".")
        tOut = tIn + ".000";
        return;
    end
    parts = split(tIn, ".");
    frac = parts(2);
    % 右侧补0到3位，超过3位截断
    if strlength(frac) < 3
        frac = frac + repmat("0", 1, 3 - strlength(frac));
    elseif strlength(frac) > 3
        frac = extractBetween(frac, 1, 3);
    end
    tOut = parts(1) + "." + frac;
end

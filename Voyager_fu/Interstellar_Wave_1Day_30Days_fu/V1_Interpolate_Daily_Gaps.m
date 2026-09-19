function [daily,audit] = V1_Interpolate_Daily_Gaps(observedDaily)
% User-authorized linear missing-day interpolation; retained for the extended figure.
% Fill only interior missing days by a straight line between real daily values.
% Explicit 'linear' in irf_resamp disables its automatic averaging branch.
% Endpoint assertions prevent its extrapolation option from ever being used.
daily=observedDaily;
observed=observedDaily.BMean_nT;
t=seconds(observedDaily.EpochUTC-observedDaily.EpochUTC(1));
missing=~isfinite(observed);
assert(all(diff(t)==86400),'Expected one original row per UTC day.');
assert(~missing(1) && ~missing(end),'Missing endpoints require extrapolation; not authorized.');
daily.BMeanObserved_nT=observed;
daily.IsInterpolated=missing;
if any(missing)
    filled=irf_resamp([t(~missing),observed(~missing)],t(missing),'linear');
    daily.BMean_nT(missing)=filled(:,2);
end
assert(isequal(daily.BMean_nT(~missing),observed(~missing)),'An observed value was changed.');
assert(all(isfinite(daily.BMean_nT)),'An internal gap was not filled.');
rows=find(missing); previous=zeros(size(rows)); following=previous;
valid=find(~missing);
for k=1:numel(rows)
    previous(k)=valid(find(valid<rows(k),1,'last'));
    following(k)=valid(find(valid>rows(k),1,'first'));
end
alpha=(t(rows)-t(previous))./(t(following)-t(previous));
expected=(1-alpha).*observed(previous)+alpha.*observed(following);
assert(all(alpha>0 & alpha<1),'Extrapolation is forbidden.');
assert(all(abs(daily.BMean_nT(rows)-expected)<1e-12),'Linear interpolation audit failed.');
audit=table(observedDaily.EpochUTC(rows),observedDaily.EpochUTC(previous), ...
    observedDaily.EpochUTC(following),observed(previous),observed(following), ...
    alpha,daily.BMean_nT(rows),following-previous-1, ...
    'VariableNames',{'InterpolatedDayUTC','PreviousObservedUTC','NextObservedUTC', ...
    'PreviousB_nT','NextB_nT','WeightOfNextValue','InterpolatedB_nT','GapLengthDays'});
end



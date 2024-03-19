function threshold =  threshold_cal(corr)
if max(corr) >= 0
    threshold = max(corr)/(10 * max(corr) - 3);
else
    threshold = max(corr) / 3;
end
end
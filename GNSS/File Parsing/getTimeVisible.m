function [periods,PRNs] = getTimeVisible(obs)
PRNs = unique(obs.PRN);
periods = cell(length(PRNs),3);

for i = 1:length(PRNs)
    ind = obs.PRN == PRNs(i);
    t = obs.tow(ind);
    dt = diff(t);
    ddt = diff(dt);
    i_i = [1; find(ddt > 0)+2];
    i_f = [find(ddt > 0)+1; length(t)];
    periods{i,1} = t(i_f) - t(i_i);
    periods{i,2} = i_i;
    periods{i,3} = i_f;
end

end

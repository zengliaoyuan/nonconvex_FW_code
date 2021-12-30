%% % plot the function value trace with respect to cpu time

%function [] = plotFvaltrace(matfilename)
clear

load('haha.mat');

for ino = 1: length(Aindex)
    maxtime = maxtime_arr(ino);
    trange = 0.01: 0.01: maxtime;
    
    errF = cell(repeat, 1);
    timeTrace = cell(repeat, 1);
    
    for rno = 1:repeat
        Fval_end = [trace_all{ino}{rno}{1}(end, 2), trace_all{ino}{rno}{2}(end, 2)];
        for kk = 1:length(Arho)
            Fval_end = [Fval_end, trace_all{ino}{rno}{2+kk}(end, 2)];
        end
        minFval = min(Fval_end);
        
        F0 = trace_all{ino}{rno}{1}(1, 2);    % note that all methods use the same initial point
        
        errF{rno} = cell(2+length(Arho), 1);
        errF{rno}{1} = (trace_all{ino}{rno}{1}(:, 2) - minFval) /  (F0 - minFval);       % errF_fw
        errF{rno}{2} = (trace_all{ino}{rno}{2}(:, 2) - minFval) /  (F0 - minFval);       % errF_afw
        for kk = 1:length(Arho)
            errF{rno}{2+kk} =  (trace_all{ino}{rno}{2+kk}(:, 2) - minFval) /  (F0 - minFval);       % err_rho{rho}
        end

        timeTrace{rno} = cell(2+length(Arho), 1);
        timeTrace{rno}{1} = trace_all{ino}{rno}{1}(:, 1);       % timetrace_fw
        timeTrace{rno}{2} = trace_all{ino}{rno}{2}(:, 1);      % timetrace_afw
        for kk = 1:length(Arho)
            timeTrace{rno}{2+kk} =  trace_all{ino}{rno}{2+kk}(:, 1);    % timetrace_rho{rho}
        end

    end
    
    FMean = cell(2+length(Arho), 1);
    
    for timeno = 1:length(trange)
        Ftall = cell(2+length(Arho), 1);
        
        for rno = 1: repeat
            Ft_fw_rno = min( errF{rno}{1}(timeTrace{rno}{1} <= trange(timeno)) );
            Ftall{1} = [Ftall{1}; Ft_fw_rno];
            
            Ft_afw_rno = min( errF{rno}{2}(timeTrace{rno}{2} <= trange(timeno)) );
            Ftall{2} = [Ftall{2}; Ft_afw_rno];
            
            for kk = 1:length(Arho)
                Ft_rho_rno = min( errF{rno}{2+kk}(timeTrace{rno}{2+kk} <= trange(timeno)) );
                Ftall{2+kk} = [Ftall{2+kk}; Ft_rho_rno];
            end
            
        end
        FMean{1} = [ FMean{1}; mean(Ftall{1}) ];            % mean E(t)_fw
        FMean{2} = [ FMean{2}; mean(Ftall{2}) ];            % mean E(t)_afw
        for kk = 1:length(Arho)
            FMean{2+kk} =  [ FMean{2+kk}; mean(Ftall{2+kk}) ] ;     % mean E(t)_rho{rho}
        end
        
    end
    
    fig = figure;
    plot(trange, FMean{1}, 'r-', 'LineWidth', 1); hold on       % plot mean E(t)_fw
    plot(trange, FMean{2}, 'm-', 'LineWidth', 1); hold on      % plot mean E(t)_afw
    color_arr = {'b',   '#D95319', 	'#7E2F8E', 'c'};
    marker_arr = {'-.' , '--'};
    for kk = 1:length(Arho)
        plot(trange, FMean{2+kk}, marker_arr{mod(kk-1, 2)+1}, 'Color', color_arr{mod(kk-1, 4)+1},  'LineWidth', 1); hold on
    end
    
    xlabel('time (s)', 'fontsize', 14);
    ylabel('Average E(t)', 'fontsize', 1)
    set(gca, 'yscale', 'log');
    titlename = ['i = ' num2str(Aindex(ino)) '  and  ' 'T^{max} = ' num2str(maxtime_arr(ino))];
    title(titlename);
    
    legend_str = cell(2+length(Arho), 1);
    legend_str(1)= {'FW'};
    legend_str(2) = {'AFW'};
    legend_str(3:end) = strcat( 'FW_{\rho_', '{', strtrim(cellstr(num2str(Arho(:)))) , '}', '}' );
    h = legend(legend_str);
    set(h, 'FontSize', 8, 'location', 'northeast')
    hold off
    
    a = clock;
    figname = ['FvalTimeTrace_index_'  num2str(Aindex(ino)) '-' date '-' num2str(a(4)) '-'  num2str(a(5))];
    savefig(fig, figname, 'Results');
end

%end

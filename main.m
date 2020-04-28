Ns = [3 7 15 31];
dts = [1/64 1/128 1/256 1/512 1/1024 1/2048 1/4096];
frames = [1/8 2/8 3/8 4/8];

for f = 1:length(frames)
    fig = figure('NumberTitle', 'off', 'Name', sprintf("Explicit: Frame = " + strtrim(rats(frames(f))) + " s"));
    for n = 1:length(Ns)
        range = 0:(1 / (Ns(n) + 1)):1;
        [X, Y] = meshgrid(range, range);
        for d = 1:length(dts)
            tic
            T = explicit_euler(Ns(n), Ns(n), dts(d), frames(f));
            toc
            T_mat = transpose(reshape(T, [(Ns(n) + 2), (Ns(n) + 2)]));
            
            
            % plot work
%             splot = subplot(length(Ns), length(dts), (length(dts) * (n - 1) + d));
%             surf(X, Y, T_mat);
%             title("Nxy=" + string(Ns(n)) + ", dt=" + strtrim(rats(dts(d))));
% 
%             seperate_fig = figure('visible','off');
%             hax_new = copyobj(splot, seperate_fig);
%             set(hax_new, 'Position', get(0, 'DefaultAxesPosition'));
%             tit = "Nxy=" + string(Ns(n)) + ", dt=" + string(dts(d))+ ", t=" + string(frames(f));
% %             saveas(seperate_fig, ('images/'+tit),'png');
%             figure(fig);
            print('done');
             
        end    
    end    
end

for f = 1:length(frames)
    fig = figure('NumberTitle', 'off', 'Name', sprintf("Implicit: Frame = " + strtrim(rats(frames(f))) + " s"));
    for n = 1:length(Ns)
        range = 0:(1 / (Ns(n) + 1)):1;
        [X, Y] = meshgrid(range, range);
        T = implicit_euler(Ns(n), Ns(n), dts(1), frames(f));
        T_mat = transpose(reshape(T, [(Ns(n) + 2), (Ns(n) + 2)]));
        
        % plot work
        subplot(2, 2, n);
        surf(X, Y, T_mat);
        title("Nxy=" + string(Ns(n)) + ", dt=" + strtrim(rats(dts(1)))); 
            
    end    
end
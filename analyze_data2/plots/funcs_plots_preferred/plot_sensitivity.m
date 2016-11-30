    
function plot_sensitivity(f,group,sens_data,chosen_sens_stat,freqband_stats,do_zscore,sr)
    
    normalize_by_sr = 0;

    if do_zscore mydata = zscore(group.data);
    else mydata = (group.data); end
    
    if normalize_by_sr
        sens_data = sens_data ./ repmat(sr(:),1,size(sens_data,2));
    end
    
    index = f >= freqband_stats(1) & f <= freqband_stats(2);
    mydata = mydata(index,:);
    mydata = mean(mydata,1);
    
    switch chosen_sens_stat
        case 1;
            sens_stat = abs(sens_data);
        case 2;
            sens_stat = [sqrt( sens_data(:,1).^2 + sens_data(:,2).^2 ) sqrt( sens_data(:,3).^2 + sens_data(:,4).^2 )];
        case 3;
            sens_stat = [ calc_ss(sens_data(:,1),sens_data(:,4)) calc_ss(sens_data(:,2),sens_data(:,3)) ];  % Note that this order is correct! See handdrawn diagram for Jefferson's ordering of ctg / Scheme's
        case 4;
            sens_stat = [ calc_ss2(sens_data) ];
        case 5
            sens_stat = [ abs(sens_data(:,1) - sens_data(:,4)) abs(sens_data(:,2) - sens_data(:,3))];
            %sens_stat(:,1) = sens_stat(:,1) ./ mean([sens_data(:,1) sens_data(:,4)],2);
            %sens_stat(:,2) = sens_stat(:,2) ./ mean([sens_data(:,2) sens_data(:,3)],2);

        case 6
            %sens_stat = [ abs(sens_data(:,1) ./ sens_data(:,4)) abs(sens_data(:,2) ./ sens_data(:,3))];
            sens_stat = log([ abs(sens_data(:,1) ./ sens_data(:,3)) abs(sens_data(:,2) ./ sens_data(:,4))]);
            sens_stat = ([ (sens_data(:,1) ./ (sens_data(:,1)+sens_data(:,3)) ) (sens_data(:,2) ./ (sens_data(:,2)+sens_data(:,4)))]);
    end
    
    %sens_stat = sqrt(sens_data(:,1).^2 + sens_data(:,2).^2);

    
    %figure; plott_fit(sens_data(:,3),mydata,'.');
    
    Nstats = size(sens_stat,2);
    Rarr = zeros(1,Nstats);
    Parr = zeros(1,Nstats);
    for i = 1:Nstats
        figure; plott_fit(sens_stat(:,i),mydata,'.');
        [R,P] = corrcoef(sens_stat(:,i),mydata);
        Rarr(i) = R(1,2); Parr(i) = P(1,2);
        xlabel('Sensitivity Statistic');
        ylabel('SFC');
    end
    figure; bar(Rarr)
    
end





function out = calc_ss(x,y)
    % Calculate the sensitivity statistic I will use
    order = 1;
    % Might need to play around with this
    out = sqrt(abs(abs(x).^order - abs(y).^order) ./ (abs(x).^order + abs(y).^order));

end



function out = calc_ss2(in)

    a=in(:,1);
    b=in(:,2);
    c=in(:,3);
    d=in(:,4);
    
    out = sqrt(abs(a.^2 + b.^2 - c.^2 - d.^2) ./ (a.^2+b.^2+c.^2+d.^2));

end
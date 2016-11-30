


function h = plot_phi_analysis(varargin)


    % Pull out inputs   
    [reg, args]=parseparams(varargin);  % Regs should contain t and X
    
    p = inputParser;
        default_plot_type = 'real_imag_scatter';
        valid_plot_types = {'real_imag_scatter','angle_cdf','angle_pdf','mag_vs_angle_scatter'};
        check_plot_type = @(x) any(validatestring(x,valid_plot_types));
    addOptional(p,'plot_type',default_plot_type,check_plot_type);
    addOptional(p,'hist_settings',{},@iscell);
    
    parse(p,args{:});
    
    % Pull out inputs
    group = reg{1};
    hist_settings = p.Results.hist_settings;
    plot_type = p.Results.plot_type;
    
    h=[];
    
    figure;
    leg_arr = {group.legend};
    for i = 1:length(group)
        switch plot_type
            case 'real_imag_scatter'

                c = get_clist(i);   % Color
                x = group(i).datastats;
                % figure; plot(abs(x),angle(x),'b.'); hold on; plot(abs(y),angle(y),'r.');
                % figure; [p,n] = ecdf(angle(x)); plot(n,p,'b'); hold on; [p,n] = ecdf(angle(y)); plot(n,p,'r');

                if i>1; figure; end     % New figure for each iteration of for loop
                hold on;
                h = plot(x,[c '.']);
                h=plot(mean(x),[c 'x'],'MarkerSize',50,'LineWidth',10);
                %viscircles([0,0],[abs(mean(x))],'EdgeColor',c);
                plotcircle(0,0,[mean(abs(x))],c);
                legend(leg_arr{i});


            case 'mag_vs_angle_scatter'
                c = get_clist(i);   % Color
                x = group(i).datastats;
                hold on;
                h=plot(abs(x),abs(angle(x)),[c '.'],'MarkerSize',20);
                ylabel('Abs(angle)'); xlabel('Mag');
                if i == length(group)
                    legend(leg_arr{:});
                end

            case 'angle_pdf'
                c = get_clist(i);   % Color
                x = group(i).datastats;
                
                if i>1; figure; end
                [n,b] = hist(abs(angle(x)),20);
                bar(b,n,c);
                xlabel('Abs(angle)'); ylabel('N');
                legend(leg_arr{i});
                
            case 'angle_cdf'
                c = get_clist(i);   % Color
                x = group(i).datastats;
                % figure; plot(abs(x),angle(x),'b.'); hold on; plot(abs(y),angle(y),'r.');
                hold on; [p,n] = ecdf(abs(angle(x))); h=plot(n,p,c);
                xlabel('Abs(angle)'); ylabel('N');
                
                if i == length(group)
                    legend(leg_arr{:});
                end
        end
    end
    
end



function plotcircle(varargin)

    xc=varargin{1};
    yc=varargin{2};
    radius=varargin{3};
    
    % Submitting linespec is optional
    if nargin < 4
        linespec = {};
    else
        linespec=varargin(4)
    end
    
    N=100;

    x = linspace(xc-radius,xc+radius,N);
    y = sqrt(radius.^2 - (x-xc).^2) + yc;
    
    hold on; plot(x,y,linespec{:});
    hold on; plot(x,-1*y,linespec{:});



end



function [ph1a, ph2a, delta_phi, Ca] = get_phasestats_ctgs(f,Cave,phiave,freqband_stats, coords, group, center180)

    if ~exist('center180','var'); center180=[]; end
    if isempty(center180); center180=1; end

    plot_on = 0;

    ind = f > freqband_stats(1) & f < freqband_stats(2);
    
    Ca = squeeze(mean(Cave(ind,:,:),1));
    ph = Cave.*exp(1i*phiave);
    ph = squeeze(mean(ph(ind,:,:),1));
    
%     ph = ph(~bad_any,:);
%     Ca = Ca(~bad_any,:);

    % Set up output coordinates - ctg0 is the absolute phase output, and
    % 1,and 2 are the differences output
    
    ctg1 = coords(1);
    ctg2 = coords(2);
    if length(coords) > 2
        ctg0 = coords(3);
    else
        ctg0 = size(Cave,3);
    end
    
    ph1 = mean(ph(:,ctg1),2);
    ph2 = mean(ph(:,ctg2),2);
    ph0 = mean(ph(:,ctg0),2);
    
    Ca1 = mean(Ca(:,ctg1),2);
    Ca2 = mean(Ca(:,ctg2),2);
    Ca = mean(Ca(:,ctg0),2);
    
    z = angle(ph1); if center180; z(z<0) = z(z<0) + 2*pi; end; ph1a = z;
    z = angle(ph2); if center180; z(z<0) = z(z<0) + 2*pi; end; ph2a = z;
    ph0a = angle(ph0); if center180; ph0a(ph0a<0) = ph0a(ph0a<0) + 2*pi; end
    
    if plot_on
        %ph0a = abs(ph0a);
        figure; hist(angle(ph0),20); xlabel('Phase (rad)'); ylabel('NCells')
        figure; hist(ph0a,20); xlabel('Phase (rad)'); ylabel('NCells')
    end
    
    
%     chosen_cells = ~bad_any(:);
%     chosen_cells = ~bad_any(:) & group(1).cells;
    chosen_cells = 1:length(ph0a);
    if plot_on
        figure; plot(ph0a(chosen_cells),Ca(chosen_cells),'k.')
        hold on; plot_scattergroups(ph0a, Ca,group(1:2));
    end
    
    
    % Difference angle
    %z = angle(ph1)-angle(ph2);
    z = ph1a-ph2a;
    ind=z>pi; z(ind)=z(ind)-2*pi;
    ind=z<-pi; z(ind)=z(ind)+2*pi;
%     %z = abs(z);
%     z = abs(ph1 - ph2);
    delta_phi = z;
    
    if plot_on

        figure; 
        hist(delta_phi,20); xlabel('ph1 - ph2');


        figure; plot(Ca,delta_phi,'k.'); xlabel('SFC');ylabel('Ph1 - Ph2'); title('Phase shifts only associated with low SFC - likely random');
        hold on; plot_scattergroups(Ca, delta_phi,group); xlabel(group(1).data_name); ylabel('Phi Diff');

        figure; plott_fit(Ca(chosen_cells),abs(delta_phi(chosen_cells)),'k.');
        hold on; plot_scattergroups(Ca, abs(delta_phi),group);
    end

end
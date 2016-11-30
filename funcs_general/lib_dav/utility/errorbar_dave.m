

function handles = errorbar_dave(t,x,e,varargin)


    handles = plot(t,x,varargin{:},'LineWidth',2); 
    hold on; plot(t,x+e,varargin{:},'LineWidth',1); 
    hold on; plot(t,x-e,varargin{:},'LineWidth',1); 

end
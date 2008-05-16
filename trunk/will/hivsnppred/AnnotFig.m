function AnnotFig(x_axis_spots,PLOT_ANNOT_CELL,FLAG)

if nargin==3
    figure
    axes

    colors=lines;
    
    num_lines=length(PLOT_ANNOT_CELL);
    
    widths=cellfun(@diff,PLOT_ANNOT_CELL(:,1));
    norm_widths=(widths./max(widths))*.8;
    
    for i=1:num_lines
        rectangle('position',[0.1 (1-(1/num_lines)*i) norm_widths(i) 1/num_lines],'facecolor',colors(rem(i,length(colors))+1,:));
        text('position',[0.1 (1-(1/num_lines)*i+0.5*(1/num_lines))],'string',PLOT_ANNOT_CELL{i,2},'fontunits','normalized','fontsize',1/num_lines)
    end
    
    set(gca,'ytick',[],'xtick',[])
    
    
    return
end

colors=lines;
bar_per=.2;

y_lim=get(gca,'ylim');
x_lim=get(gca,'xlim');
set(gca,'ylim',[y_lim(1)-y_lim(2)*bar_per y_lim(2)])

if x_lim(1)>=0
    width=x_lim(2)-x_lim(1);
    x_lim=[x_lim(1)-0.05*width x_lim(2)+0.05*width];
    set(gca,'xlim',x_lim)
end



total_bar_height=bar_per*y_lim(2);
ind_bar_height=total_bar_height/5;

binary_mask=false(size(PLOT_ANNOT_CELL,1),length(x_axis_spots));

for k=1:size(PLOT_ANNOT_CELL,1)
    spots=sort(PLOT_ANNOT_CELL{k,1});
    first_spot=find(x_axis_spots>spots(1),1);
    second_spot=find(x_axis_spots>spots(2),1);
    if isempty(second_spot)
        trans_spots=[length(x_axis_spots) x_lim(2)];
        
        
    elseif first_spot==1&&second_spot==1
        trans_spots=[x_lim(1) 0];
    elseif second_spot==first_spot
        trans_spots=[first_spot second_spot+1];
    else
        %trans_spots=[x_axis_spots(first_spot) x_axis_spots(second_spot)];
        trans_spots=[first_spot second_spot];
        
    end

    plot_row=find(sum(binary_mask(:,first_spot:second_spot),2)==0,1);
    binary_mask(plot_row,first_spot:second_spot)=true;

    rectangle('position',[trans_spots(1) -plot_row*ind_bar_height trans_spots(2)-trans_spots(1) ind_bar_height],'facecolor',colors(rem(k,length(colors))+1,:));
    text('position',[trans_spots(1)+(trans_spots(2)-trans_spots(1))/2 -plot_row*ind_bar_height+.5*ind_bar_height],'string',PLOT_ANNOT_CELL{k,2},'fontunits','normalized','fontsize',.01)

end

end
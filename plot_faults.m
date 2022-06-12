function [sfl]=plot_faults(faults,dim,xlim,ylim,col)
%function [sfl]=plot_faults(faults,dim,xlim,ylim,col);
%
%fid=fopen('ECSZ.flt','wt');
[sfl,l]=size(dim);
for i=1:sfl-1
    if i==1
        S1=1;
    else
        S1=S2+1;
    end
    S2=S1+dim(i)-1;
    ab=faults(S1:S2,1);
    or=faults(S1:S2,2);
    if min(ab)>xlim(1) & max(ab)<xlim(2) & min(or)>ylim(1) & max(or)<ylim(2)
        %    line(ab,or,'Color',col,'LineStyle','-','LineWidth',1.5,'HandleVisibility','off'), hold on
        line(ab,or,'Color',col,'LineStyle','-','LineWidth',1.0), hold on
        %    plot3(ab,or,0*ab,'Color',col,'LineStyle','-','LineWidth',2,'HandleVisibility','off'), hold on
        %    line(ab,or,'Color',[0.2 0.2 0.2],'LineStyle','-','LineWidth',1), hold on
        %    line(ab,or,'Color','w','LineStyle','-','LineWidth',1), hold on
        %    fprintf(fid,'> \n');
        %    for j=1:length(ab)
        %      fprintf(fid,'%9.4f %8.4f \n',ab(j),or(j));
        %    end
    end
end
%fclose(fid);


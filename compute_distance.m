function d = compute_distance()

clear d
load sa_meg_template
para.grid = 'medium';

pos = tp_grid2aal(sa_meg_template.grid_medium',para);

for i = 1 : 90
  for j = 1 : 90
    
    d(i,j) = ((pos(1,i)-pos(1,j))^2+(pos(2,i)-pos(2,j))^2+(pos(3,i)-pos(3,j))^2)^1/2;
    
  end
end

d(d==inf)=0;

% figure_white; plot3(pos(1,:),pos(2,:),pos(3,:),'k.')
% 
% for i = 1 : 90
%   for j = 1 : 90
%     
%     linewidth = floor(SC(i,j)*30);
%     if ~linewidth < 1
%     line([pos(1,i) pos(1,j)],[pos(2,i) pos(2,j)],[pos(3,i) pos(3,j)],'linewidth',linewidth,'color',[rand rand rand])
%     end
%   end
% end
%     
%     
%     
    








% s=trisurf(sa.vc_indi.tri,sa.vc_indi.vc(:,1),sa.vc_indi.vc(:,2),sa.vc_indi.vc(:,3),0.5*ones([size(sa.vc_indi.vc(:,1),1) 1]))
% hold on
% h=light; material([0.4 0.6 0.2]); lighting GOURAUD
% shading interp;

% rotate(s,[0 0 1],270)

% plot3(sa.locs_3D_indi(:,1),sa.locs_3D_indi(:,2),sa.locs_3D_indi(:,3),'.')

% 
% x = sa.grid_fine_indi(:,1);
% y = sa.grid_fine_indi(:,2);
% z = sa.grid_fine_indi(:,3);
% 
% xv = linspace(min(x),max(x),100);
% yv = linspace(min(y),max(y),100);
% 
% [X,Y]=meshgrid(xv,yv);
% 
% Z = griddata(x,y,z,X,Y);

close all
s=trisurf(sa.segmentation{3}.tri,sa.segmentation{3}.vc(:,1),sa.segmentation{3}.vc(:,2),sa.segmentation{3}.vc(:,3),0.5*ones([size(sa.segmentation{3}.vc(:,1),1) 1]))
 hold on
h=light; material([0.5 0.5 0.5]); lighting GOURAUD
shading interp;
alpha(0.3); axis off

s=trisurf(sa.cortex10K.tri,sa.cortex10K.vc(:,1),sa.cortex10K.vc(:,2),sa.cortex10K.vc(:,3),0.5*ones([size(sa.cortex10K.vc(:,1),1) 1]))
h=light; material([0.5 0.5 0.5]); lighting GOURAUD

shading interp;


a = sa.locs_3D_indi*sa.trafo.u_indi2template+sa.trafo.r_indi2template
plot3(a(:,1),a(:,2),a(:,3),'.','color',[0.6 0.6 0.6],'markersize',12)














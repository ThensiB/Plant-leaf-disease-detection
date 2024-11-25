function colorTransformationStructure(Cmd)

figure(2)
switch Cmd
case 1
    subplot(1,2,1)
    [x,y,z,Tri] = makeshape('Cube');
    CData = [x,y,z];
    myplot((x-0.5)*0.8,(y-0.5)*0.8,(z)*0.8,Tri,CData);
    coloraxis('R',5,0.5*0.8);
    coloraxis('G',6,0.5*0.8); 
    coloraxis('B',3);
case 2
    subplot(1,2,2)
    [x,y,z,Tri] = makeshape('Double Hexacone');
    CData = colorspace('rgb<-hsi',[(pi+atan2(-y,-x))*180/pi,sqrt(x.^2+y.^2),z]);      
    myplot(x,y,2*z,Tri,CData);
    coloraxis('H',1);
    coloraxis('S',2);   
    coloraxis('I',4);
end

axis equal;
axis off;
pbaspect([1,1,1]);
view(70,27);
rotate3d on;

return;


function myplot(x,y,z,Tri,CData)
% Plot a triangular mesh with color data
cla;
CData = min(max(CData,0),1);
patch('Faces',Tri,'Vertices',[x,y,z],'FaceVertexCData',CData,...
   'FaceColor','interp','EdgeColor','none');
hold on;
return;


function coloraxis(Name,Type,h)
% Draw color axes as 3D objects
FontSize = 14;

switch Type
case 1
   set(text(-0.25,-1.3,1.1,Name),'FontWeight','bold',...
      'FontSize',FontSize,'Color',[0,0,0]);
   t = linspace(0,pi*3/2,60);
   x = cos(t)*1.1;
   y = sin(t)*1.1;
   set(plot3(x,y,zeros(size(x)),'k-'),'LineWidth',2.5,'Color',[1,1,1]*0.8);
   x = [x,-0.1,0,-0.1];
   y = [y,-1.05,-1.1,-1.15];
   set(plot3(x,y,ones(size(x)),'k-'),'LineWidth',2.5);   
case 2
   set(text(0,-0.6,0.15,Name),'FontWeight','bold','FontSize',FontSize,'Color',[0,0,0]);
   set(plot3([-0.05,0,0.05,0,0],[-0.9,-1,-0.9,-1,0],[0,0,0,0,0],'k-'),'LineWidth',2.5);
case 3
   set(text(0,0.15,1.3,Name),'FontWeight','bold','FontSize',FontSize,'Color',[0,0,0]);
   set(plot3([0,0,0,0,0],[-0.05,0,0.05,0,0],[1.6,1.7,1.6,1.7,0],'k-'),'LineWidth',2.5);
case 4
   set(text(0,0.15,2.4,Name),'FontWeight','bold','FontSize',FontSize,'Color',[0,0,0]);
   set(plot3([0,0,0,0,0],[-0.05,0,0.05,0,0],[2.4,2.5,2.4,2.5,0],'k-'),'LineWidth',2.5);
case 5
   set(text(0.6,0,h+0.15,Name),'FontWeight','bold','FontSize',FontSize,'Color',[0,0,0]);
   set(plot3([0.9,1,0.9,1,-1,-0.9,-1,-0.9],[-0.05,0,0.05,0,0,-0.05,0,0.05],...
      zeros(1,8)+h,'k-'),'LineWidth',2.5);
case 6
   set(text(0,0.6,h+0.15,Name),'FontWeight','bold','FontSize',FontSize,'Color',[0,0,0]);
   set(plot3([-0.05,0,0.05,0,0,-0.05,0,0.05],[0.9,1,0.9,1,-1,-0.9,-1,-0.9],...
      zeros(1,8)+h,'k-'),'LineWidth',2.5);
end
return;


function [x,y,z,Tri] = makeshape(Shape)

N = 12;    
Nth = 25;
Nz = 8;  

switch lower(Shape)
case 'cube'
   [u,v] = meshgrid(linspace(0,1,N),linspace(0,1,N));
   u = u(:);
   v = v(:);
   x = [u;u;u;u;zeros(N^2,1);ones(N^2,1)];
   y = [v;v;zeros(N^2,1);ones(N^2,1);v;v];
   z = [zeros(N^2,1);ones(N^2,1);v;v;u;u];
   Tri = trigrid(N,N);
   Tri = [Tri;N^2+Tri;2*N^2+Tri;3*N^2+Tri;4*N^2+Tri;5*N^2+Tri];
case 'double hexacone'
   Nz = floor(Nz/2)*2+1;
   [u,v] = meshgrid(linspace(0,2*pi,Nth),linspace(0,1,Nz));
   Tri = trigrid(Nth,Nz);
   r = 1 - abs(2*v - 1);
   r = 0.92*r./max(max(abs(cos(u)),abs(cos(u - pi/3))),abs(cos(u + pi/3)));
   x = r(:).*cos(u(:));
   y = r(:).*sin(u(:));
   z = v(:);   
end
return;


function Tri = trigrid(Nu,Nv)
i = (1:(Nu-1)*Nv-1).';
i(Nv:Nv:end) = [];
Tri = [i,i+1,i+Nv;i+1,i+Nv+1,i+Nv];
return;
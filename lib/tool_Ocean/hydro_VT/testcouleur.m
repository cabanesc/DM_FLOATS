% Teste les couleurs des palettes

colormap(jet);

z=ones(65,1)*[0:64];
v=0:65;

figure
contourf(z,v)
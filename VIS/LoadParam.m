% Load all the simulation paraemters
compe=double(h5read('param','/compe'));
resgrid=double(h5read('param','/resgrid'));
c=double(h5read('param','/c'));
mi=double(h5read('param','/mi'));
me=double(h5read('param','/me'));
Nflvr=h5read('param','/Nflvr');
epc=double(h5read('param','/epc'));
nx=h5read('param','/nx');
ny=h5read('param','/ny');
nz=h5read('param','/nz');
psaveratio=h5read('param','/psaveratio');
extBx=double(h5read('param','/extBx'));
extBy=double(h5read('param','/extBy'));
extBz=double(h5read('param','/extBz'));
extBmag=sqrt(extBx.^2+extBy.^2+extBz.^2);

sx=double(nx/resgrid)+1;
sy=double(ny/resgrid)+1;
sz=double(nz/resgrid)+1;
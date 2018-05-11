%% this program provides a quick visualisation of a field varaible in 3D
DataFolder='/Users/rahul/Dropbox/NumericalCodes/PICTOR/data';
FldName='/Jx2';
Xslice=1;
Yslice=100;
Zslice=10;
tfld=24000:24000:24000;%Default time steps to load field data
%% =======================================================================
CurrentFolder=pwd;
cd(DataFolder)
compe=double(h5read('param','/compe'));
resgrid=double(h5read('param','/resgrid'));
c=double(h5read('param','/c'));

i=1;
for t=tfld
FileName=strcat('fld_',num2str(t)); 
fld=h5read(FileName,FldName);
%fld=smooth3(fld,'box',[2 2 2]);

%fldx=h5read(FileName,'/Bx'); 
%fldy=h5read(FileName,'/By'); 
%fldz=h5read(FileName,'/Bz'); 
%fld=fldx.^2+fldy.^2+fldz.^2; 


%plot the field intensity in the slices 
h=slice(fld,Xslice,Yslice,Zslice);
colorbar
set(h,'edgeColor','none')
title(strcat('Time =',strcat(num2str(floor(t*c/compe)),'\omega_{pe}^{-1}')));
xlabel('c/\omega_{pe}')
ylabel('c/\omega_{pe}')
zlabel('c/\omega_{pe}')
axis equal
pause(0.0)
drawnow;

end 

cd(CurrentFolder)
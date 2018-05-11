%% a quick visualisation of the field qunatities
DataFolder='/Users/rahulkumar/Dropbox/NumericalCodes/PICTOR/data';
TimeSteps=100:100:4800;
FldName='/Bz';
%% ==============================================================================
CurrentFolder=pwd;
cd(DataFolder)

compe=double(h5read('param','/compe'));
resgrid=double(h5read('param','/resgrid'));
c=double(h5read('param','/c'));

for t=TimeSteps
    t
    str=t*c/compe;
    FileName=strcat('fld_',num2str(t));
    fld=h5read(FileName,FldName);
    %fld=sign(fld).*abs(fld).^(0.25);
    if(~exist('h1'))
        fig=figure('units','normalized','outerposition',[0 0 1 1]) ;
    end
        h1=surf(transpose(fld),'EdgeColor','none','facecolor','interp');
        grid off
        axis equal
        
        xlabel('2x / (c\omega_{pe}^{-1})','FontSize',14);
        ylabel('2y / (c\omega_{pe}^{-1})','FontSize',14);
        
        colorbar
        view([0 0 1])
    
        caxis([-max(abs(fld(:))) max(abs(fld(:)))]);
        title(strcat('t=',strcat(num2str(str),'\omega_{pe}^{-1}')))
        drawnow;
end


cd(CurrentFolder);
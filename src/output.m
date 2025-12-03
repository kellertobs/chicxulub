if ~mod(step,nout)
    if lvplt
        VIS = {'Visible','on'};
    else
        VIS = {'Visible','off'};
    end

    if ~exist('fh2','var'); fh2 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh2); clf;
    end

    colormap(ocean);

    axh = 6.00; axw = 7.50; %   Height and width of axis
    ahs = 1.50; avs = 1.20; %   Horzontal and vertial distance between axis
    axb = 1.50; axt = 2.00; %   Bottom and top;Size of page relative to axis
    axl = 1.75; axr = 0.90; %   Right and left; spacing of axis to page

    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh2,UN{:},'Position',[9 9 fw fh]);
    set(fh2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh2,'Color','w','InvertHardcopy','off');
    set(fh2,'Resize','off');
    ax(1) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(2) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(3) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
    ax(4) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(5) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    ax(6) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

    sgtitle(sprintf('Time elapsed %.3f yr', time/3600/24/365.25),TX{:},FS{:})

    set(fh2, 'CurrentAxes', ax(1))
    imagesc(x,z,-w.*3600); axis equal tight; box on; cb = colorbar; hold on;
    if ~any(isnan(unit(:)))
        for i = 1:size(unit,3)
            contour(x,z,unit(:,:,i),1,'w','LineWidth',0.5);
        end
    end
    clim([min(-w(:).*3600),max(-w(:).*3600)])
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Segregation z-speed [m/hr]',TX{:},FS{:})
    ylabel('Depth [m]',TX{:},FS{:});

    set(fh2, 'CurrentAxes', ax(2))
    imagesc(x,z,u.*3600); axis equal tight;  box on; cb = colorbar; hold on;
    if ~any(isnan(unit(:)))
        for i = 1:size(unit,3)
            contour(x,z,unit(:,:,i),1,'w','LineWidth',0.5);
        end
    end
    clim([min(u(:).*3600),max(u(:).*3600)])
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Segregation x-speed [m/hr]',TX{:},FS{:})

    set(fh2, 'CurrentAxes', ax(3))
    imagesc(x,z,p(2:end-1,2:end-1)); axis equal tight;  box on; cb = colorbar; hold on;
    if ~any(isnan(unit(:)))
        for i = 1:size(unit,3)
            contour(x,z,unit(:,:,i),1,'w','LineWidth',0.5);
        end
    end
    clim([min(p(:)),max(p(:))])
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Dynamic fluid pressure [Pa]',TX{:},FS{:})

    set(fh2, 'CurrentAxes', ax(4))
    imagesc(x,z,T); axis equal tight; box on; cb = colorbar; hold on
    if ~any(isnan(unit(:)))
        for i = 1:size(unit,3)
            contour(x,z,unit(:,:,i),1,'w','LineWidth',0.5);
        end
    end
    clim([min(Tin(:)-eps),max(Tin(:)+eps)])
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Temperature [$^\circ$C]',TX{:},FS{:})
    ylabel('Depth [m]',TX{:},FS{:});
    xlabel('Distance [m]',TX{:},FS{:});

    set(fh2, 'CurrentAxes', ax(5))
    imagesc(x,z,C); axis equal tight;  box on; cb = colorbar; hold on;
    if ~any(isnan(unit(:)))
        for i = 1:size(unit,3)
            contour(x,z,unit(:,:,i),1,'w','LineWidth',0.5);
        end
    end
    clim([min(Cin(:)-eps),max(Cin(:)+eps)])
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Salinity [wt]',TX{:},FS{:})
    xlabel('Distance [m]',TX{:},FS{:});

    set(fh2, 'CurrentAxes', ax(6))
    imagesc(x,z,V); axis equal tight;  box on; cb = colorbar; hold on
    if ~any(isnan(unit(:)))
        for i = 1:size(unit,3)
            contour(x,z,unit(:,:,i),1,'w','LineWidth',0.5);
        end
    end
    clim([min(V(:)-eps),max(V(:)+eps)])
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Vapour [wt]',TX{:},FS{:})
    xlabel('Distance [m]',TX{:},FS{:});
    drawnow

    % print figure and save data to file
    if svout
        print(fh2,[outdir,'/',runID,'/',runID,'_sol_',int2str(step/nout)],'-dpng','-r200')
        save([outdir,'/',runID,'/',runID,'_',int2str(step/nout)],'x','z','u','w','p','f','T','C','V','dTdt','dCdt','dVdt','K','Drho','time','Ra','unit');
    end
end
function [F] = DHEM_plot(FWR_times,FWR,Atlas,X,Y,Z,in,slices,argin)

nz = size(X,1);
nx = size(X,2);
ny = size(X,3);

ncell = nz*ny*nx;

y = reshape(Y(1,1,:),ny,1);
interp_scheme= argin;
% Plot interp fields
set(figure(100), 'Position', [25 100 1800 900])

nslices = length(slices);
nsbplot = nslices * 6;
    
for tt = 1 : length(FWR_times)

    % Interp Forward Field    
    FWR_dBx_interp = griddata(FWR.XYZ(:,3),FWR.XYZ(:,1),FWR.XYZ(:,2),FWR.dBx(:,tt),...
        reshape(Z,ncell,1),reshape(X,ncell,1),reshape(Y,ncell,1),interp_scheme);
    
    % FWR_dBx_interp(isnan(FWR_dBx_interp)==1) = 0;
    % save([work_dir '\' 'dBdt_interp\dBxdt\dBx_t' num2str(tt) '.dat'],'-ascii','dBx_dt_interp');
    
    FWR_dBx_interp = reshape(FWR_dBx_interp,nz,nx,ny);
%     FWR_dBx_interp(FWR_dBx_interp>0) = FWR_dBx_interp(FWR_dBx_interp>0)/max(FWR_dBx_interp(FWR_dBx_interp>0));
%     FWR_dBx_interp(FWR_dBx_interp<0) = FWR_dBx_interp(FWR_dBx_interp<0)/abs(min(FWR_dBx_interp(FWR_dBx_interp<0)));
    
    % Plot interp fields    
    for kk = 1:nslices
        subplot(nslices,6,(nsbplot-1) - 6*(kk-1)-3);
        imagesc(FWR_dBx_interp(:,:,slices(kk)));
        caxis([min(min(min(FWR_dBx_interp))) max(max(max(FWR_dBx_interp)))]);axis off
        caption = ['\bfY:' num2str(y(slices(kk)))];
        text(1,4,caption,'Color',[1 1 1]); hold on;
        if kk == nslices
            title('\bfSQUID dBx/dt');
        end

        if kk == 1
            axis on
            caption = ['\bfTime:' num2str(FWR_times(tt))];
            xlabel(caption);
        end

    end
        
    FWR_dBy_interp = griddata(FWR.XYZ(:,3),FWR.XYZ(:,1),FWR.XYZ(:,2),FWR.dBy(:,tt),...
    reshape(Z,ncell,1),reshape(X,ncell,1),reshape(Y,ncell,1),interp_scheme);
    
    % FWR_dBy_interp(isnan(FWR_dBy_interp)==1) = 0;
    % save([work_dir '\' 'dBdt_interp\dBxdt\dBx_t' num2str(tt) '.dat'],'-ascii','dBx_dt_interp');
    
    FWR_dBy_interp = reshape(FWR_dBy_interp,nz,nx,ny);
%     FWR_dBy_interp(FWR_dBy_interp>0) = FWR_dBy_interp(FWR_dBy_interp>0)/max(FWR_dBy_interp(FWR_dBy_interp>0));
%     FWR_dBy_interp(FWR_dBy_interp<0) = FWR_dBy_interp(FWR_dBy_interp<0)/abs(min(FWR_dBy_interp(FWR_dBy_interp<0)));
    
    for kk = 1:nslices
        subplot(nslices,6,(nsbplot-2) - 6*(kk-1));
        imagesc(FWR_dBy_interp(:,:,slices(kk)));
        caxis([min(min(min(FWR_dBy_interp))) max(max(max(FWR_dBy_interp)))]);axis off
        caption = ['\bfY:' num2str(y(slices(kk)))];
        text(1,4,caption,'Color',[1 1 1]); hold on;
        if kk == nslices
            title('\bfSQUID dBy/dt');
        end

        if kk == 1
            axis on
            caption = ['\bfTime:' num2str(FWR_times(tt))];
            xlabel(caption);
        end

    end
        
    FWR_dBz_interp = griddata(FWR.XYZ(:,3),FWR.XYZ(:,1),FWR.XYZ(:,2),FWR.dBz(:,tt),...
    reshape(Z,ncell,1),reshape(X,ncell,1),reshape(Y,ncell,1),interp_scheme);
    
    % FWR_dBz_interp(isnan(FWR_dBz_interp)==1) = 0;
    % save([work_dir '\' 'dBdt_interp\dBxdt\dBx_t' num2str(tt) '.dat'],'-ascii','dBx_dt_interp');
    
    FWR_dBz_interp = reshape(FWR_dBz_interp,nz,nx,ny);
%     FWR_dBz_interp(FWR_dBz_interp>0) = FWR_dBz_interp(FWR_dBz_interp>0)/max(FWR_dBz_interp(FWR_dBx_interp>0));
%     FWR_dBz_interp(FWR_dBz_interp<0) = FWR_dBz_interp(FWR_dBz_interp<0)/abs(min(FWR_dBz_interp(FWR_dBz_interp<0)));
    
    for kk = 1:nslices
        subplot(nslices,6,(nsbplot) - 6*(kk-1));
        imagesc(FWR_dBz_interp(:,:,slices(kk)));
        caxis([min(min(min(FWR_dBz_interp))) max(max(max(FWR_dBz_interp)))]);axis off
        caption = ['\bfY:' num2str(y(slices(kk)))];
        text(1,4,caption,'Color',[1 1 1]); hold on;
        if kk == nslices
            title('\bfSQUID dBz/dt');
        end

        if kk == 1
            axis on
            caption = ['\bfTime:' num2str(FWR_times(tt))];
            xlabel(caption);
        end

    end
        
    % dBx Select specific time channels and interp
%     index = abs(t_channel-FWR_times(tt));
%     index = find(index==min(index));
    % index =tt;
    dBx_interp = griddata(Atlas.DXYZ(in,4),Atlas.DXYZ(in,2),Atlas.DXYZ(in,3),Atlas.dBx(in,tt),...
        reshape(Z,ncell,1),reshape(X,ncell,1),reshape(Y,ncell,1),interp_scheme);

    % dBx_dt_interp(isnan(dBx_dt_interp)==1) = 0;


    dBx_interp = reshape(dBx_interp,nz,nx,ny);
%     dBx_interp(dBx_interp>0) = dBx_interp(dBx_interp>0)/max(dBx_interp(dBx_interp>0));
%     dBx_interp(dBx_interp<0) = dBx_interp(dBx_interp<0)/abs(min(dBx_interp(dBx_interp<0)));


    for kk = 1:nslices
        subplot(nslices,6,(nsbplot-2) - 6*(kk-1)-3);
        imagesc(dBx_interp(:,:,slices(kk)));
        caxis([min(min(min(FWR_dBx_interp))) max(max(max(FWR_dBx_interp)))]);axis off
        caption = ['\bfY:' num2str(y(slices(kk)))];
        text(1,4,caption,'Color',[1 1 1]); hold on;
        if kk == nslices
            title('\bfDATA dBx/dt');
        end

        if kk == 1
            axis on
            caption = ['\bfTime:' num2str(FWR_times(tt))];
            xlabel(caption);
        end

    end

    dBy_interp = griddata(Atlas.DXYZ(in,4),Atlas.DXYZ(in,2),Atlas.DXYZ(in,3),Atlas.dBy(in,tt),...
        reshape(Z,ncell,1),reshape(X,ncell,1),reshape(Y,ncell,1),interp_scheme);

    % dBy_dt_interp(isnan(dBy_dt_interp)==1) = 0;


    dBy_interp = reshape(dBy_interp,nz,nx,ny);
%     dBy_interp(dBy_interp>0) = dBy_interp(dBy_interp>0)/max(dBy_interp(dBy_interp>0));
%     dBy_interp(dBy_interp<0) = dBy_interp(dBy_interp<0)/abs(min(dBy_interp(dBy_interp<0)));

    for kk = 1:nslices
        subplot(nslices,6,(nsbplot) - 6*(kk-1)-3);
        imagesc(dBy_interp(:,:,slices(kk)));
        caxis([min(min(min(FWR_dBy_interp))) max(max(max(FWR_dBy_interp)))]);axis off
        caption = ['\bfY:' num2str(y(slices(kk)))];
        text(1,4,caption,'Color',[1 1 1]); hold on;
        if kk == nslices
            title('\bfDATA dBy/dt');
        end

        if kk == 1
            axis on
            caption = ['\bfTime:' num2str(FWR_times(tt))];
            xlabel(caption);
        end

    end

    dBz_interp = griddata(Atlas.DXYZ(in,4),Atlas.DXYZ(in,2),Atlas.DXYZ(in,3),Atlas.dBz(in,tt),...
        reshape(Z,ncell,1),reshape(X,ncell,1),reshape(Y,ncell,1),interp_scheme);

    % dBz_dt_interp(isnan(dBz_dt_interp)==1) = 0;


    dBz_interp = reshape(dBz_interp,nz,nx,ny);
%     dBz_interp(dBz_interp>0) = dBz_interp(dBz_interp>0)/max(dBz_interp(dBz_interp>0));
%     dBz_interp(dBz_interp<0) = dBz_interp(dBz_interp<0)/abs(min(dBz_interp(dBz_interp<0)));

    for kk = 1:nslices
        subplot(nslices,6,(nsbplot-1) - 6*(kk-1));
        imagesc(dBz_interp(:,:,slices(kk)));
        caxis([min(min(min(FWR_dBz_interp))) max(max(max(FWR_dBz_interp)))]);axis off
        caption = ['\bfY:' num2str(y(slices(kk)))];
        text(1,4,caption,'Color',[1 1 1]); hold on;
        if kk == nslices
            title('\bfDATA dBz/dt');
        end

        if kk == 1
            axis on
            caption = ['\bfTime:' num2str(FWR_times(tt))];
            xlabel(caption);
        end

    end


        
    F(tt) = getframe(gcf);
    fprintf('Time: %i done!\n',tt);
    figure(100);hold off

end





 

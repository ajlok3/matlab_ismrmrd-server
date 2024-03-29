classdef perfusion_rest < handle
    methods
        function process(obj, connection, config, metadata, logging)
            logging.info('Config: \n%s', config);
            zf_counter = 0;
            rep_counter = 0;
            % Metadata should be MRD formatted header, but may be a string
            % if it failed conversion earlier
            try
                logging.info("Incoming dataset contains %d encodings", numel(metadata.encoding))
                logging.info("First encoding is of type '%s', with field of view of (%g x %g x %g)mm^3, matrix size of (%g x %g x %g), and %g coils", ...
                    metadata.encoding(1).trajectory, ...
                    metadata.encoding(1).encodedSpace.fieldOfView_mm.x, ...
                    metadata.encoding(1).encodedSpace.fieldOfView_mm.y, ...
                    metadata.encoding(1).encodedSpace.fieldOfView_mm.z, ...
                    metadata.encoding(1).encodedSpace.matrixSize.x, ...
                    metadata.encoding(1).encodedSpace.matrixSize.y, ...
                    metadata.encoding(1).encodedSpace.matrixSize.z, ...
                    metadata.acquisitionSystemInformation.receiverChannels)
            catch
                logging.info("Improperly formatted metadata: \n%s", metadata)
            end

            % Continuously parse incoming data parsed from MRD messages
            acqGroup = cell(1,0); % ismrmrd.Acquisition;
            imgGroup = cell(1,0); % ismrmrd.Image;
            wavGroup = cell(1,0); % ismrmrd.Waveform;
            try
                while true
                    item = next(connection);

                    % ----------------------------------------------------------
                    % Raw k-space data messages
                    % ----------------------------------------------------------
                    if isa(item, 'ismrmrd.Acquisition')
                        % Accumulate all imaging readouts in a group
                        if (~item.head.flagIsSet(item.head.FLAGS.ACQ_IS_NOISE_MEASUREMENT)    && ...
                            ~item.head.flagIsSet(item.head.FLAGS.ACQ_IS_PHASECORR_DATA)       && ...
                            ~item.head.flagIsSet(item.head.FLAGS.ACQ_IS_PARALLEL_CALIBRATION)       )
                                acqGroup{end+1} = item;
                        end
                        
                        % When this criteria is met, run process_raw_zero_filled() on the accumulated
                        % data, which returns images that are sent back to the client.
                        if (item.head.flagIsSet(item.head.FLAGS.ACQ_LAST_IN_REPETITION))
                            rep_counter = rep_counter + 1;
                            if mod(rep_counter, 2) == 0
                                group_length = length(acqGroup) - zf_counter;
                                zf_counter = length(acqGroup);
                                logging.info(sprintf("Processing a group of k-space data; %d items", length(acqGroup(end-group_length+1:end))))
                                zero_filled = obj.process_raw_zero_filled(acqGroup(end-group_length+1:end), config, metadata, logging);
                                logging.debug("Sending zero-filled images to client")
                                connection.send_image(zero_filled);
                            end
                        end

                    elseif isempty(item)
                        break;

                    else
                        logging.error("Unhandled data type: %s", class(item))
                    end
                end
            catch ME
                logging.error(sprintf('%s\nError in %s (%s) (line %d)', ME.message, ME.stack(1).('name'), ME.stack(1).('file'), ME.stack(1).('line')));
            end

            
            % Process group of raw k-space data.
            if ~isempty(acqGroup)
                logging.info("Processing a group of k-space data (untriggered)")
                 image = obj.process_raw(acqGroup, config, metadata, logging);
                 logging.debug("Sending image to client")
                 connection.send_image(image);
                
                acqGroup = cell(1,0);
            end

            connection.send_close();
            return
        end
        
        % Process a set of raw k-space data and return the zero-filled image
        function images = process_raw_zero_filled(obj, group, config, metadata, logging)
            images = cell(1,0);

            % This is almost like the twix_obj
            twix_obj = twix_map_obj_fire;
            twix_obj.setMrdAcq(group);
            
            %% this is necessary to avoid non-square images and 'oversized' kspace
            twix_obj.NLin = twix_obj.NCol / 2;
            twix_obj.NRep = 1;
            twix_obj.Rep = ones(1,length(twix_obj.Rep));
            
            %% preparing kspace
            kspAll = squeeze(twix_obj.imageData());
            if length(size(kspAll)) > 3
                kspAll = kspAll(:,:,:,2); % take second slice, if more than one available
            end
            kspAll = permute(kspAll, [1,3,2]);
            
            %% oversampling reduction
            tmpKsp =fftshift(fft(fftshift(kspAll,1),[],1),1)/sqrt(size(kspAll,1));
            tmpKsp = tmpKsp(size(tmpKsp, 1)/4 + 1 : size(tmpKsp, 1)*3/4, :, :, :, :, :, :, :, :);
            kspAll = fftshift(ifft(fftshift(tmpKsp,1),[],1),1)*sqrt(size(tmpKsp,1));

           %% sum of squares iFFT recon
            ref=ifft2c_mri(kspAll);
            ref = sum(abs(ref),3);  
            
            % image processing
            % TODO: export complex images
            img = abs(ref);
            logging.debug("Image data is size %d x %d", size(img))

            % Normalize and convert to short (int16)
            img = img .* (32767./max(img(:)));
            img = int16(round(img));
            img = rot90(img, 2);
            % Invert image contrast
            % img = int16(abs(32767-img));

            % Format as ISMRMRD image data
            % TODO: send back slice by slice
            for sl=1:size(img, 4)
                for tm=1:size(img, 3)
                    image = ismrmrd.Image(img(:, :, tm, sl));

                    % Find the center k-space index
                    % TODO: fix center Idx
                    centerIdx = find((twix_obj.Lin == twix_obj.centerLin) & (twix_obj.Sli == sl), 1);
                    if isempty(centerIdx)
                        warning('Could not find center k-space readout')
                        centerIdx = 1;
                    end
                    % centerIdx = 1;

                    % Copy the relevant AcquisitionHeader fields to ImageHeader
                    image.head = image.head.fromAcqHead(group{centerIdx}.head);
                    
                    
                    % field_of_view is mandatory
                    image.head.field_of_view  = single([metadata.encoding(1).reconSpace.fieldOfView_mm.x ...
                                                        metadata.encoding(1).reconSpace.fieldOfView_mm.y ...
                                                        metadata.encoding(1).reconSpace.fieldOfView_mm.z]);
                    
                    % calib data ==> series index 0                                
                    image.head.image_series_index = 0;
                    image.head.image_index = tm;
                    
                    % Set ISMRMRD Meta Attributes
                    meta = struct;
                    meta.DataRole               = 'Image';
                    meta.ImageProcessingHistory = 'MATLAB';
                    meta.WindowCenter           = uint16(16384);
                    meta.WindowWidth            = uint16(32768);
                    meta.ImageRowDir            = group{centerIdx}.head.read_dir;
                    meta.ImageColumnDir         = group{centerIdx}.head.phase_dir;
                    meta.InstanceNumber         = tm;
                    meta.ImageComments          = 'ZF_recon';
                    meta.SiemensControl_SkipSaveOnHost = {'bool', 'true'};
                                        
                    % set_attribute_string also updates attribute_string_len
                    image = image.set_attribute_string(ismrmrd.Meta.serialize(meta));
                    images{end+1} = image;
                end
            end
            logging.info(sprintf('Reconstructed %d zero-filled images', numel(images)))
        end


        % Process a set of raw k-space data and return an image
        function images = process_raw(obj, group, config, metadata, logging)
            images = cell(1,0);

            % This is almost like the twix_obj
            twix_obj = twix_map_obj_fire;
            twix_obj.setMrdAcq(group);
            
            kspAll = twix_obj.imageData();
            tmpKsp =fftshift(fft(fftshift(kspAll,1),[],1),1)/sqrt(size(kspAll,1));
            tmpKsp = tmpKsp(size(tmpKsp, 1)/4 + 1 : size(tmpKsp, 1)*3/4, :, :, :, :, :, :, :, :);
            kspAll = fftshift(ifft(fftshift(tmpKsp,1),[],1),1)*sqrt(size(tmpKsp,1));
            logging.info("Data is 'mapVBVD formatted' with dimensions:")  % Data is 'mapVBVD formatted' with dimensions:
            logging.info(sprintf('%s ', twix_obj.dataDims{1:10}))         % Col Cha Lin Par Sli Ave Phs Eco Rep Set
            logging.info(sprintf('%3d ', size(kspAll)))                   % 404  14 124   1   1   1   1   1   1  11

            %%
            %  Daming's code from here
            %%
            kspAll = squeeze(kspAll);
            if length(size(kspAll)) < 5
                kspTmp(:,:,:,1,:) = kspAll;
                kspAll = kspTmp;
                clear kspTmp
            end
            kspAll = permute(kspAll, [1, 3, 5, 2, 4]);
            lambda1 = 0.010;
            N_PD = 9;
            
            [Nx, Ny, Nt, Nc, Nslc] = size(kspAll);
            for slices = 1:Nslc
                % slices = 1
                coil_keep = 0;
                Recon_CS_TTV_TPCA = [];
                kSpace_slc = kspAll(:,:,:,:,slices);

                % detecting noisy coils
                CS_mask=zeros(size(kSpace_slc,2),size(kSpace_slc,3));
                for zz=1:size(CS_mask,2)
                    CS_mask(find(kSpace_slc(size(kSpace_slc,1)/2+1,:,zz,1)),zz)=1;
                end
                
                CS_mask2 = zeros(Nx,Ny,Nt);
                
                for i = 1:Nx
                    CS_mask2(i,:,:) = CS_mask;
                end
                
                count = 1;
                for coil = 1:Nc
                    kSpace_tmp = squeeze(kspAll(:,:,1,coil,1)).*squeeze(CS_mask2(:,:,1));
                    [col,row,kSpace_new] = find(abs(kSpace_tmp));
                    [N,edges] = histcounts(kSpace_new,88);
                    maxCounts = max(N);
                    leftBin = find(N > maxCounts/2, 1, 'first');
                    rightBin = find(N > maxCounts/2, 1, 'last');
                    fwhm = rightBin - leftBin;
                    if fwhm < 10
                        coil_keep(count) = coil;
                        count = count +1;
                    end
                 end



                kSpace_slc = kspAll(:,:,:,coil_keep, slices);
                [Nx, Ny, Nt, Nc] = size(kSpace_slc);

                 % % Coil Compression: PCA 8
                disp('Peforming PCA on coils...');
                NumOfComp = 8;
                 if Nc > NumOfComp
                    tmp_data = double( squeeze(kSpace_slc(:,:,:,:)) );
                    tmp_data = reshape( tmp_data, [Nx*Ny*Nt Nc] );
                    [coeff, score, latent] = pca( tmp_data );
                    compressed_data = tmp_data*coeff(:,1:NumOfComp);
                    compressed_data = reshape( compressed_data, [Nx Ny Nt NumOfComp] );
                    kSpace_slc = compressed_data;
        %             Coil_kdata_ref = kSpace_slc;
                 else
                     NumOfComp = Nc;
                 end
                    % clear tmp_data compressed_data coeff



                    % move PD to the end and scaling PD scans
                    kSpace_slc_temp = kSpace_slc(:,:,[N_PD+1:Nt,1:N_PD],:);
                    clear kSpace_slc;
                    kSpace_slc = kSpace_slc_temp;
                    t1w = squeeze(kSpace_slc(:,:,Nt-N_PD,:));
                    pd = squeeze(kSpace_slc(:,:,Nt-N_PD+1,:));
                    tmp_sort_t1w = sort( abs(t1w(:)), 'descend' );
                    tmp_sort_pd = sort( abs(pd(:)), 'descend' );
                    mean_t1w = mean(tmp_sort_t1w(1:round(length(t1w(:))*1)));
                    mean_pd = mean(tmp_sort_pd(1:round(length(pd(:))*1)));
                    max_ratio = mean_pd/mean_t1w
                    kSpace_slc(:,:,Nt-N_PD+1:Nt,:) = kSpace_slc(:,:,Nt-N_PD+1:Nt,:)./max_ratio;
                    % move PD to the end and scaling PD scans

                    CS_mask=zeros(size(kSpace_slc,2),size(kSpace_slc,3));
                    for zz=1:size(CS_mask,2)
                        CS_mask(find(kSpace_slc(size(kSpace_slc,1)/2+1,:,zz,1)),zz)=1;
                    end
                    
                    CS_mask2 = zeros(Nx,Ny,Nt);
                    for i = 1:Nx
                        CS_mask2(i,:,:) = CS_mask;
                    end
        % % %             DC=squeeze(sum(kSpace_slc,3));
                    DC = squeeze( sum(kSpace_slc(:,:,1:end,:),3) );
                    SM=squeeze(sum(CS_mask2,3));
                    SM(find(SM==0))=1;
                    for ii = 1:NumOfComp
                        DC(:,:,ii)=DC(:,:,ii)./SM;
                    end
                    %%% This ref is the unpaired FFT for the coil sensitivity maps. The
                    %%% unpaired FFT for the kdata is withinEmat_2DPERF.mtimes
                    ref=ifft2c_mri(DC);
                    [dummy,b1]=adapt_array_2d_st2(ref);
                    b1=b1/max( abs(b1(:)) );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % iterative reconstruction
                    param = struct();
                    param.TPCA = TempPCA();
                    param.TTV = TV_Temp();
                    param.STV = TV_2DPERF();
                    param.nite = 9;
                    param.display = 1;
                    param.VW = 'on';
                    param.TTVw = lambda1;
                    param.TPCAw = lambda1*0.5;
                    %param.TPCAw = 0;
                    param.STVw = 0;

                    param.E = Emat_2DPERF( gpuArray(single(CS_mask2)), gpuArray(single(b1)) );
                    %recon the first image
                    param.y = gpuArray(single( kSpace_slc ));
                    Recon_CS_TTV_TPCA = gpuArray( param.E'*param.y );

                    param.block_size = 8;
                    param.tau_size = 0.15;

                    ITER = 3; %3
                    for ite = 1:ITER
                        disp( [ sprintf( 'iteration = %2d', ite )] )
                        Recon_CS_TTV_TPCA = CSL1NlCg_TTV_TPCA_STV_VW_gpuArray( Recon_CS_TTV_TPCA, param );
                    end
                    Recon_CS_TTV_TPCA_CoilCompGPU = gather( Recon_CS_TTV_TPCA );
                    % % scaling BACK
                    Recon_CS_TTV_TPCA_CoilCompGPU(:,:,Nt-N_PD+1:Nt) = Recon_CS_TTV_TPCA_CoilCompGPU(:,:,Nt-N_PD+1:Nt).*max_ratio;
                    recon_cs_total(:,:,:,slices) = Recon_CS_TTV_TPCA_CoilCompGPU(:,:,[Nt-N_PD+1:Nt,1:Nt-N_PD]);
            end
            % image processing
            % TODO: export complex images
            img = abs(recon_cs_total);
            logging.debug("Image data is size %d x %d x %d after coil combine and phase oversampling removal", size(img))

            % Normalize and convert to short (int16)
            img = img .* (32767./max(img(:)));
            img = int16(round(img));
            img = rot90(img, 2);

            % Format as ISMRMRD image data
            % TODO: send back slice by slice
            %arr = {};
            for sl=1:size(img, 4)
                for tm=1:size(img, 3)
                    image = ismrmrd.Image(img(:, :, tm, sl));

                    % Find the center k-space index
                    % TODO: fix center Idx
                    centerIdx = find((twix_obj.Lin == twix_obj.centerLin) & (twix_obj.Sli == sl) & (twix_obj.Rep == tm), 1);
                    if isempty(centerIdx)
                        warning('Could not find center k-space readout')
                        centerIdx = 1;
                    end
                    % centerIdx = 1;

                    % Copy the relevant AcquisitionHeader fields to ImageHeader
                    image.head = image.head.fromAcqHead(group{centerIdx}.head);
                    
%                     mil_sec = image.head.acquisition_time_stamp;
%                     res = double(mod(mil_sec, 400))/400;
%                     sec = floor(double(mil_sec)/400);
%                     mins = floor(sec/60);
%                     hours = floor(mins/60);
%                     res = res + hours*10000 + mod(mins,60)*100 +mod(sec,60);
%                     disp(res);
%                     arr = [arr; [res, sl, tm]]; 
                    % field_of_view is mandatory
                    image.head.field_of_view  = single([metadata.encoding(1).reconSpace.fieldOfView_mm.x ...
                                                        metadata.encoding(1).reconSpace.fieldOfView_mm.y ...
                                                        metadata.encoding(1).reconSpace.fieldOfView_mm.z]);
                    
                                                    
                    image.head.image_series_index = sl;
                    image.head.image_index = tm;
                    
                    % Set ISMRMRD Meta Attributes
                    meta = struct;
                    meta.DataRole               = 'Image';
                    meta.ImageProcessingHistory = 'MATLAB';
                    meta.WindowCenter           = uint16(16384);
                    meta.WindowWidth            = uint16(32768);
                    meta.ImageRowDir            = group{centerIdx}.head.read_dir;
                    meta.ImageColumnDir         = group{centerIdx}.head.phase_dir;
                    meta.InstanceNumber         = tm;
                    meta.ImageComments = metadata.measurementInformation.protocolName;
                    if tm <= N_PD
                        meta.SiemensDicom_adFlipAngleDegree = {'double', '5.0'};
                    end
                    if sl == 1
                        meta.SequenceDescriptionAdditional      = '_AIF';
                        meta.ImageComments = strcat(meta.ImageComments, '_AIF');
                    end
                   
                    
                    
                    % set_attribute_string also updates attribute_string_len
                    image = image.set_attribute_string(ismrmrd.Meta.serialize(meta));

                    images{end+1} = image;
                end
            end
%             dlmwrite('test.csv', arr, 'delimiter', ',', 'precision', 12);
            logging.info(sprintf('Reconstructed %d images', numel(images)))
        end

        % Placeholder function that returns images without modification
        function images = process_images(obj, group, config, metadata, logging)
            images = group;
        end
    end
end

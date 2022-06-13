classdef cine_flash_cartesian < handle
    methods
        function process(obj, connection, config, metadata, logging)
            logging.info('Config: \n%s', config);
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
                        
                        % discard the first 300 measurement before steady
                        % state is reached
                        
                        % Accumulate all imaging readouts in a group
                        if (~item.head.flagIsSet(item.head.FLAGS.ACQ_IS_NOISE_MEASUREMENT)    && ...
                            ~item.head.flagIsSet(item.head.FLAGS.ACQ_IS_PHASECORR_DATA)       && ...
                            ~item.head.flagIsSet(item.head.FLAGS.ACQ_IS_PARALLEL_CALIBRATION)       )
                                acqGroup{end+1} = item;
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

        % Process a set of raw k-space data and return an image
        function images = process_raw(obj, group, config, metadata, logging)
            images = cell(1,0);

            % This is almost like the twix_obj
            twix_obj = twix_map_obj_fire;
            twix_obj.setMrdAcq(group);
            
            kspAll = twix_obj.imageData();
            % crop first 25 timeframes before steady state is reached
            kspAll = kspAll(:,:,:,:,:,:,26:end);
            tmpKsp =fftshift(fft(fftshift(kspAll,1),[],1),1)/sqrt(size(kspAll,1));
            tmpKsp = tmpKsp(size(tmpKsp, 1)/4 + 1 : size(tmpKsp, 1)*3/4, :, :, :, :, :, :, :, :);
            kspAll = fftshift(ifft(fftshift(tmpKsp,1),[],1),1)*sqrt(size(tmpKsp,1));
            logging.info("Data is 'mapVBVD formatted' with dimensions:")  % Data is 'mapVBVD formatted' with dimensions:
            logging.info(sprintf('%s ', twix_obj.dataDims{1:10}))         % Col Cha Lin Par Sli Ave Phs Eco Rep Set
            logging.info(sprintf('%3d ', size(kspAll)))                   % 404  14 124   1   1   1   1   1   1  11

            kspAll = squeeze(kspAll);
            if length(size(kspAll)) < 5
                kspTmp(:,:,:,1,:) = kspAll;
                kspAll = kspTmp;
                clear kspTmp
            end
            kspAll = permute(kspAll, [1, 3, 5, 2, 4]);
            
            lambda1 = 0.010;
            N_PD = 9;
            NumOfCoilCmp = 8;
            
            [Nsample,Nlin,Nphs,Ncoil,Nslc] = size(kspAll);
            for SLC = 1:Nslc
                fprintf( 'Slice = %02d\n', SLC )        

                kSpace_SLC = kspAll(:,:,:,:,SLC);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % CS mask
                CS_mask = abs(squeeze(kSpace_SLC(:,:,:,1)));
                CS_mask( find( CS_mask > 0 ) ) = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Excluding noise coils (low SNR)
        if 1
            DC = squeeze(sum(kSpace_SLC,3));
            SM = squeeze(sum(CS_mask,3));
            SM(find(SM==0)) = 1;
            for Nc = 1 : size(kSpace_SLC,4)
                DC(:,:,Nc) = DC(:,:,Nc)./SM;
            end
            ref = ifft2c_mri(DC);
            ref = ref./max(abs(ref(:)));
            CoilInfo = [];
            for Nc = 1:size(ref,3)
                tmp_input = abs( ref(:,:,Nc) );
                CoilInfo(1,Nc) = mean( tmp_input(:) );
                CoilInfo(2,Nc) = std( tmp_input(:) ); 
            end
            SignalRatio = CoilInfo(2,:)./CoilInfo(1,:);
            ind = find( SignalRatio < 1.0 );
            kSpace_SLC( :,:,:,ind ) = [];
        end
        
        Ncoil2 = size(kSpace_SLC,4);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coil compression: PCA to make 8 virtual coils
        NumOfCoils = NumOfCoilCmp;
        if Ncoil2 > NumOfCoils
            disp('Peforming PCA for coil compression...');
            data = reshape(kSpace_SLC,[Nsample*Nlin*Nphs  Ncoil2]);
            clear kSpace_SLC
            [coeff,~,~] = pca(data); %principal component analysis
            compressed_data = data*coeff(:,1:NumOfCoils);
            kSpace_SLC = reshape(compressed_data,[Nsample Nlin Nphs NumOfCoils]);
        else
            disp('No coil compression...');
            NumOfCoils = Ncoil2;
        end
        clear data compressed_data coeff
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coil sensitivity map
disp('Coil sens map...');
        if 0
            DC = squeeze(sum(kSpace_SLC,3));
            SM = squeeze(sum(CS_mask,3));
        else
            StartCoilSensPhase = 6;
            DC = squeeze(sum(kSpace_SLC(:,:,StartCoilSensPhase:end,:),3));
            SM = squeeze(sum(CS_mask(:,:,StartCoilSensPhase:end,:),3));
        end
        SM(find(SM==0)) = 1;
        for ii = 1 : NumOfCoils
            DC(:,:,ii) = DC(:,:,ii)./SM;
        end
        %%% This ref is the unpaired FFT for the coil sensitivity maps.
        %%% The unpaired FFT for the kdata is within Emat_2DPERF.mtimes
        ref = ifft2c_mri(DC);
        [~,CoilSens] = adapt_array_2d(ref);
        CoilSens = CoilSens ./ max(abs(CoilSens(:)));
        clear  DC  SM 
        clear  ref

        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CS - preparation
        clear param
        param.TPCA = TempPCA();
        param.TTV = TV_Temp();
        param.STV = TV_2DPERF();
        param.nite = 8;
        param.display = 1;
        param.STVw = 0;
        param.Slice = SLC;

        % % GPU mode
        param.E = Emat_2DPERF( single(CS_mask), single(CoilSens) );
        param.y = single(kSpace_SLC);
        
        ITER = 4;
        param.VW = 'on';

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CS - iterative recon
        lambda =  0.01; 
        param.TTVw = lambda;
        param.TPCAw = param.TTVw/10;
        Recon_IterCS = gpuArray( param.E'*param.y );              

        tic
        for ite = 1 : ITER
            fprintf( 'CS Rrecon: Slice = %02d,  Iteration = %d out of %d\n', SLC, ite,  ITER )
            % % GPU mode
            Recon_IterCS = CSL1NlCg_TTV_TPCA_STV_VW_gpuArray( Recon_IterCS, param );
            toc
        end
        CSrecon_all(:,:,:,SLC) = gather(Recon_IterCS);  
        
    end  % slice
            % image processing
            % TODO: export complex images
            img = abs(CSrecon_all);
            logging.debug("Image data is size %d x %d x %d after coil combine and phase oversampling removal", size(img))

            % Normalize and convert to short (int16)
            img = img .* (32767./max(img(:)));
            img = int16(round(img));
            img = rot90(img, 2);

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
                    
                    
                    % set_attribute_string also updates attribute_string_len
                    image = image.set_attribute_string(ismrmrd.Meta.serialize(meta));

                    images{end+1} = image;
                end
            end
            logging.info(sprintf('Reconstructed %d images', numel(images)))
        end

        % Placeholder function that returns images without modification
        function images = process_images(obj, group, config, metadata, logging)
            images = group;
        end
    end
end

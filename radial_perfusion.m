classdef radial_perfusion < handle
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
            %% parameter setting
            AdvDCFb1 = 1;
            
            images = cell(1,0);

            % This is almost like the twix_obj
            twix_obj = twix_map_obj_fire;
            twix_obj.setMrdAcq(group);
            
            % this is necessary to avoid non-square images and 'oversized' kspace
            twix_obj.NRep = 2;
            twix_obj.Rep = twix_obj.Rep - min(twix_obj.Rep) + 1;
           
            kspAll = twix_obj.imageData();
%             tmpKsp =fftshift(fft(fftshift(kspAll,1),[],1),1)/sqrt(size(kspAll,1));
%             tmpKsp = tmpKsp(size(tmpKsp, 1)/4 + 1 : size(tmpKsp, 1)*3/4, :, :, :, :, :, :, :, :);
%             kspAll = fftshift(ifft(fftshift(tmpKsp,1),[],1),1)*sqrt(size(tmpKsp,1));
            logging.info("Data is 'mapVBVD formatted' with dimensions:")  % Data is 'mapVBVD formatted' with dimensions:
            logging.info(sprintf('%s ', twix_obj.dataDims{1:10}))         % Col Cha Lin Par Sli Ave Phs Eco Rep Set
            logging.info(sprintf('%3d ', size(kspAll)))                   % 404  14 124   1   1   1   1   1   1  11

            
            %% sorting dimensions in correct manner.
            kSpace = permute(squeeze(kspAll),[1 3 4 5 2]);
            kSpace = kSpace(:,:,:,end-1:end,:);
            [nx,ny,nz,nt,nc] = size(kSpace);
            clear angle_info
            for i = 1:length(group)
               angle_info(i) = double(group{i}.head.user_int(5));
            end
            angle_info = reshape(angle_info, [ny, nz, nt])./10000;
            angle_3 = squeeze(angle_info(:,1,:));
            

            [nx,ny,nz,nt,nc] = size(kSpace);

            iz = min(nz,2);

            
            kSpace_SLC = squeeze(kSpace(:,:,iz,:,:));
            clear tmp_kSpace  tmp_angleinfo  ray_index  recon_nufft 
            ray_index = 1:ny;
            tmp_angleinfo = [];
            tmp_kSpace = [];
            tmp_kSpace = reshape( kSpace_SLC, [nx, ny*nt, nc] );
            tmp_angleinfo = col( angle_3 );

            ray_grid = -0.5:1/nx:0.5-1/nx;
            ray_info = [];
            for aaa = 1:nt*ny
                ray_info(:,aaa) = ray_grid.*exp(1i*tmp_angleinfo(aaa));
            end
            ku = 1i.*(-real(ray_info)) - imag(ray_info);
            sel_DCF = 2;   % 1 = standard DCF % 2 = new DCF
            switch sel_DCF
                case 1
                     wu = abs( ku )./max(abs(ku(:)));
                case 2
                    Option.Nsamples = nx;
                    Option.AnglePrecision = 4;   % full precision
                    Option.Display = 0;
                    Option.WeightedContrast = 0;
                    if AdvDCFb1 == 1
                       wu = AdvancedDCF_2DRadial( tmp_angleinfo, Option );
                    end
            end

            param.E = MCNUFFT_2D_GPU_kxyt_single(ku, wu, ones(nx,nx,nc));       
            param.y = tmp_kSpace;
            NUFFT_CoilSens = squeeze( param.E'*param.y);
            ref = sum(abs(NUFFT_CoilSens),3);
            
            % image processing
            % TODO: export complex images
            img = abs(ref);
            logging.debug("Image data is size %d x %d", size(img))

            % Normalize and convert to short (int16)
            img = img .* (32767./max(img(:)));
            img = int16(round(img));
            img = rot90(img, 2)';
            % Invert image contrast
            % img = int16(abs(32767-img));

            % Format as ISMRMRD image data
            % TODO: send back slice by slice
            image = ismrmrd.Image(img);

            % Find the center k-space index
            centerIdx = find((twix_obj.Lin == twix_obj.centerLin) & (twix_obj.Sli == iz), 1);
            if isempty(centerIdx)
                warning('Could not find center k-space readout')
                centerIdx = 1;
            end

            % Copy the relevant AcquisitionHeader fields to ImageHeader
            image.head = image.head.fromAcqHead(group{centerIdx}.head);


            % field_of_view is mandatory
            image.head.field_of_view  = single([metadata.encoding(1).reconSpace.fieldOfView_mm.x ...
                                                metadata.encoding(1).reconSpace.fieldOfView_mm.y ...
                                                metadata.encoding(1).reconSpace.fieldOfView_mm.z]);

            % calib data ==> series index 0                                
            image.head.image_series_index = 0;
            image.head.image_index = 0;

            % Set ISMRMRD Meta Attributes
            meta = struct;
            meta.DataRole               = 'Image';
            meta.ImageProcessingHistory = 'MATLAB';
            meta.WindowCenter           = uint16(16384);
            meta.WindowWidth            = uint16(32768);
            meta.ImageRowDir            = group{centerIdx}.head.read_dir;
            meta.ImageColumnDir         = group{centerIdx}.head.phase_dir;
            meta.InstanceNumber         = 0;
            meta.ImageComments          = 'ZF_recon';
            meta.SiemensControl_SkipSaveOnHost = {'bool', 'true'};

            % set_attribute_string also updates attribute_string_len
            image = image.set_attribute_string(ismrmrd.Meta.serialize(meta));
            images{end+1} = image;
            logging.info(sprintf('Reconstructed %d zero-filled images', numel(images)))
        end


        % Process a set of raw k-space data and return an image
        function images = process_raw(obj, group, config, metadata, logging)
            %% parameter setting

            % coil compression: if we want to do coil compression, set to 1
            ccFlag = 1;

            % If want to use pyo's DCF to calculate b1 sensitivity map, set to 1
            AdvDCFb1 = 1;
            
            images = cell(1,0);

            % This is almost like the twix_obj
            twix_obj = twix_map_obj_fire;
            twix_obj.setMrdAcq(group);
            
            kspAll = twix_obj.imageData();
%             tmpKsp =fftshift(fft(fftshift(kspAll,1),[],1),1)/sqrt(size(kspAll,1));
%             tmpKsp = tmpKsp(size(tmpKsp, 1)/4 + 1 : size(tmpKsp, 1)*3/4, :, :, :, :, :, :, :, :);
%             kspAll = fftshift(ifft(fftshift(tmpKsp,1),[],1),1)*sqrt(size(tmpKsp,1));
            logging.info("Data is 'mapVBVD formatted' with dimensions:")  % Data is 'mapVBVD formatted' with dimensions:
            logging.info(sprintf('%s ', twix_obj.dataDims{1:10}))         % Col Cha Lin Par Sli Ave Phs Eco Rep Set
            logging.info(sprintf('%3d ', size(kspAll)))                   % 404  14 124   1   1   1   1   1   1  11

            
            %% sorting dimensions in correct manner.
            kSpace = permute(squeeze(kspAll),[1 3 4 5 2]); 
            [nx,ny,nz,nt,nc] = size(kSpace);
            
            for i = 1:length(group)
               angle_info(i) = double(group{i}.head.user_int(5));
            end
            angle_info = reshape(angle_info, [ny, nz, nt])./10000;
            angle_3 = squeeze(angle_info(:,1,:));
            
            %% trajectory corrections for both FP and PD
            % kSpace: input kSpace data with (nx,nray,ncoil),
            %         nray should be time-combined
            % param.Indexes: kSpace ray number that to use for gradient delay estimation,
            %                better pick rays that are at steady state: e.g., Nray-100:Nray
            % param.N : angle scheme; e.g., tiny golden angle N=5
            % param.AngleRange :  1) '-H' for [0,180)
            %                     2) '-G' for [0,360)
            kSpace_TF = kSpace(:,1:end,:,1:end-2,:);
            [nx,ny,nz,nt,nc] = size(kSpace_TF);
            angle_3_TF = angle_3(1:end,1:end-2);
            
            for iz = 1:nz
            clear tmp_kSpace kSpace_SLC tmp_angle
            kSpace_SLC = squeeze(kSpace_TF(:,:,iz,:,:));
            tmp_kSpace = reshape( kSpace_SLC, [nx, ny*nt, nc] );
            tmp_angle = col(angle_3_TF);

            Option_SLC.Indexes = (ny*nt+1-600):ny*nt;  % last 600 rays
            Option_SLC.N = 5;
            Option_SLC.AngleRange = '-G';
            Option_SLC.Angles = tmp_angle;
            if isfield( Option_SLC, 'PC_coeff' )
               Option_SLC = rmfield( Option_SLC, 'PC_coeff' );
            end
            disp( 'Trajectory correction (RING)...' )
            Option_SLC.PC_coeff = RING_TrajCorrCoefficient( tmp_kSpace, Option_SLC );
            Option_SLC.PhaseCorrection = 'UseTC';
            PCcorrected = RadialTrajCorrection_WithRING( kSpace_SLC, kSpace_SLC, Option_SLC );
            kSpace_SLC = PCcorrected.kSpace;
            clear PCcorrected
            kSpace_tc_TF(:,:,iz,:,:) = kSpace_SLC;
            end    

            % kSpace_new(:,15:42,:,:,:) = kSpace_tc_TF;
            kSpace_new = kSpace_tc_TF;
            angle_3_new = angle_3;
            %%
            clear kSpace angle_3
            % clear kSpace
            [nx,ny,nz,nt,nc] = size(kSpace_new);
            kSpace = kSpace_new;%% slcie 4/5/6
            %kSpace = kSpace_new(:,1:ny,:,:,[2:10,13,15,16,19,21,23,24,26:30]); %% slcie 1/2/3


            angle_3 = angle_3_new(1:ny,1:end-2);

            [nx,ny,nz,nt,nc] = size(kSpace);
            if ccFlag == 1
               % PCA to make 8 virtual coils
               disp('Peforming PCA on coils...');
               no_comp = 8;
               data = reshape(kSpace,[nx*ny*nz*nt nc]);
               clear kSpace
               [coeff,~,~] = pca(data); %principal component analysis
               compressed_data = data*coeff(:,1:no_comp);
               kSpace = reshape(compressed_data,[nx ny nz nt no_comp]);
               clear data compressed_data coeff
            end

            [nx,ny,nz,nt,nc] = size(kSpace);

            %% b1 sensitivity map
            b1 = complex(zeros(nx,nx,nz,nc));

            for iz = 1:nz
                kSpace_SLC = squeeze(kSpace(:,:,iz,:,:));
                clear tmp_kSpace  tmp_angleinfo  ray_index  recon_nufft 
                ray_index = 1:ny;
                tmp_angleinfo = [];
                tmp_kSpace = [];
                tmp_kSpace = reshape( kSpace_SLC, [nx, ny*nt, nc] );
                tmp_angleinfo = col( angle_3 );

                fprintf( 'Calculate coil sensitivity....\n' );
                if iz == 1
                ray_grid = -0.5:1/nx:0.5-1/nx;
                ray_info = [];
                for aaa = 1:nt*ny
                    ray_info(:,aaa) = ray_grid.*exp(1i*tmp_angleinfo(aaa));
                end
                ku = 1i.*(-real(ray_info)) - imag(ray_info);
                sel_DCF = 2;   % 1 = standard DCF % 2 = new DCF
                switch sel_DCF
                    case 1
                         wu = abs( ku )./max(abs(ku(:)));
                    case 2
                        Option.Nsamples = nx;
                        Option.AnglePrecision = 4;   % full precision
                        Option.Display = 0;
                        Option.WeightedContrast = 0;
                        if AdvDCFb1 == 1
                           wu = AdvancedDCF_2DRadial( tmp_angleinfo, Option );
                        end
                end
            %         wu = ones(nx, ny*nt);
                end

                param.E = MCNUFFT_2D_GPU_kxyt_single( ku, wu, ones(nx,nx,nc) );       
                param.y = tmp_kSpace;
                NUFFT_CoilSens = squeeze( param.E'*param.y );
                %figure(102); imagescn( abs( NUFFT_CoilSens ), [0 0.00005],[],[],3 )
                [dummy,b1sens] = adapt_array_2d_st2( NUFFT_CoilSens );
                b1sens = b1sens./max( abs(b1sens(:)) );
                % figure(104); imagescn( abs( dummy), [],[],[],3 )
                tmp_nufft(:,:,iz) = sqrt( sum( NUFFT_CoilSens.^2, 3 ) );
                % figure(101); imagescn( abs( b1sens ),[],[],[] );   
                b1(:,:,iz,:) = b1sens;
                NUFFT_CoilSens_all(:,:,iz,:) = NUFFT_CoilSens;
                dummy_all(:,:,iz) = dummy;
                clear b1sens NUFFT_CoilSens dummy
            end

            %%
            Ray_keep = 30;
            [nx, ny, nz, nt, nc]=size(kSpace);
            clear kSpace_new angle_4
            kSpace_new = kSpace(:,ny-Ray_keep+1:ny,:,:,:);
            angle_4 = angle_3(ny-Ray_keep+1:ny,:);
            clear kSpace angle_3
            kSpace = kSpace_new;
            angle_3 = angle_4;
            [nx, ny, nz, nt, nc]=size(kSpace);
            %%
            half_width = 3; % half_width_test(wid);
            WeightedContrast = ones( nx, ny, nt);
            hourglass_filt_mask_use = WeightedContrast;

            %% CS recon
            for iz = 1:nz
                kSpace_SLC = squeeze(kSpace(:,:,iz,:,:));
                clear tmp_kSpace  tmp_angleinfo  ray_index  recon_nufft
                ray_index = 1:ny;
                tmp_angleinfo = [];
                tmp_kSpace = [];
                tmp_kSpace = kSpace_SLC;
                tmp_angleinfo = ( angle_3 );
                ray_grid = -0.5:1/nx:0.5-1/nx;
                ray_info = [];
                for aaa = 1:nt
                     for bbb = 1:ny
                         ray_info(:,bbb,aaa) = ray_grid.*exp(1i*tmp_angleinfo(bbb,aaa));
                     end
                end
                ku = 1i.*(-real(ray_info)) - imag(ray_info);

                hourglass_filt_mask_use (find(hourglass_filt_mask_use == 0)) = eps;
                w_temp = squeeze(hourglass_filt_mask_use(:,:,1));
                wu = hourglass_filt_mask_use;


                b1sens = squeeze(b1(:,:,iz,:));
                param.E = MCNUFFT_2D_GPU_kxyt( ku, wu, b1sens );       
                param.y = permute(tmp_kSpace, [1,2,4,3]);
                Recon_RadPerf_grid = param.E'*param.y;

                Weight = 0.005;
                param.TV = TV_Temp();
                times_TV = 1;
                param.TVw = Weight*times_TV;%.*max(abs(dummy1(:)));
                times_PCA = 0;
                param.PCA = TempPCA();   
                param.PCAw = param.TVw .* times_PCA;
                times_STV = 0;
                param.L1w = Weight.*times_STV;
                param.W = TV_2DPERF();

                param.nite = 9;
                param.display = 1;
            %     param.mask = mask_temp;

                recon_cs = squeeze(Recon_RadPerf_grid);
                ite = 3;
                for n = 1 : ite
                    recon_cs = CSL1NlCg_TV_PCA_new_TA_test_VW(recon_cs,param);
                end   
                recon_all(:,:,:,iz) = recon_cs;
            end
            % image processing
            % TODO: export complex images
            img = abs(recon_all);
            logging.debug("Image data is size %d x %d x %d after coil combine and phase oversampling removal", size(img))

            % Normalize and convert to short (int16)
            img = img .* (32767./max(img(:)));
            img = int16(round(img));
            img = permute(rot90(img,2), [2, 1, 3:length(size(img))]);

            % Format as ISMRMRD image data
            % TODO: send back slice by slice
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
                    if sl == 1
                        meta.SequenceDescriptionAdditional      = '_AIF';
                        meta.ImageComments = strcat(meta.ImageComments, '_AIF');
                    end
                   
                    
                    
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

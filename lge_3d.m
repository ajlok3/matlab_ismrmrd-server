classdef lge_3d < handle
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
                            if 0 % mod(rep_counter, 2) == 0
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
            %% recon params
            sel_gpu =  1;         % 1 or 2
            gpuDevice( sel_gpu )
            SEL_ImgNavi     = 0;  %ImgNavigator: 0=off; 1=on
            SEL_VS          = 1;  %ViewSharing:  0=off; 1=on
            SEL_SaveData    = 0;  %Save data: 0=off; 1=on
            DisplayFigure   = 'off';  % on or off
    
            NaviSetting = [6 4 1 1   10.0  10.0   15 ];
            SEL_PCA = 2;  
            Index_ExcludingCoils =  [ ]     
            RotatingAngle = -90;
            Index_ExcludingCoils_Navi = [  ]
            CoilTreatment = 0;
            images = cell(1,0);
            
            % This is almost like the twix_obj
            twix_obj = twix_map_obj_fire;
            twix_obj.setMrdAcq(group);
            
            kspAll = twix_obj.imageData();
            kspAll = permute(kspAll, [1,3,4,2]);
            [Nasymmetric, Nshot, Npartition, Ncoil] = size(kspAll);
            AsymmetryCenter = group{2}.head.center_sample;
            angle_mdb = twix_obj.iceParam(5,:)';
            angle_mdb = angle_mdb./10000;
            PartInfo = twix_obj.iceParam(6,:)';

            logging.info("Data is 'mapVBVD formatted' with dimensions:")  % Data is 'mapVBVD formatted' with dimensions:
            logging.info(sprintf('%s ', twix_obj.dataDims{1:10}))         % Col Cha Lin Par Sli Ave Phs Eco Rep Set
            logging.info(sprintf('%3d ', size(kspAll)))                   % 404  14 124   1   1   1   1   1   1  11

            %%
            %  Pyo's code from here
            %%
            Nsample = (Nasymmetric - AsymmetryCenter)*2;
            MaxScanInfo = max(PartInfo(:));  

            kSpace = zeros( Nsample, Nshot, Npartition, Ncoil);
            AngleInfo = zeros(Nshot, Npartition );
            
            Navigator = []; Angle_Navi = [];

            for iM = 1:length(group)
                % % below is for the self-navi at the beginning
                lin_iM = group{iM}.head.idx.kspace_encode_step_1 + 1;
                par_iM = group{iM}.head.idx.kspace_encode_step_2 + 1;
                raw_iM = group{iM}.data;

                if PartInfo(iM) <= 1  
                    Navigator(Nsample-Nasymmetric+1:Nsample,lin_iM, PartInfo(iM)+1,:) = raw_iM;
                    Angle_Navi( lin_iM, PartInfo(iM)+1 ) = angle_mdb(iM);
                elseif  PartInfo(iM) >=  MaxScanInfo-5 
                    Navigator(Nsample-Nasymmetric+1:Nsample,lin_iM, PartInfo(iM)-19, : ) = raw_iM;
                    Angle_Navi( lin_iM, PartInfo(iM)-19 ) = angle_mdb(iM);
                else
                    kSpace(Nsample-Nasymmetric+1:Nsample,lin_iM,par_iM, : ) = raw_iM;
                    AngleInfo(lin_iM,par_iM) = (angle_mdb(iM));
                end
                if   PartInfo(iM) == (MaxScanInfo-4)/2  
                    Navigator(Nsample-Nasymmetric+1:Nsample, lin_iM, PartInfo(iM)-9, : ) = raw_iM;
                    Angle_Navi( lin_iM, PartInfo(iM)-9 ) = angle_mdb(iM);
                end
            end
            clear nMeasurements  lin  par  phs subLine  cha   respiration
            clear rawdata2  iM    angle_mdb  mdh
        
        [Nsample, Nshot, Npartition,  Ncoil] = size( kSpace );
        size(kSpace)
        size(Navigator)
        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% % correction for delayed trajectories by using RING method
% % kSpace: input kSpace data with (nx,nray,ncoil),
% %         nray should be time-combined
% % param.Indexes: kSpace ray number that to use for gradient delay estimation,
% %                better pick rays that are at steady state: e.g., Nray-100:Nray
% % param.N : tiny golden angle N (e.g., 5).
% % param.AngleRange :  1) '-H' for [0,180)
% %                     2) '-G' for [0,360)
fprintf( 'Trajectory correction in each coil...\n' )    
        kSpace_TC = kSpace;
            CorrectionRange = 2:Npartition;
            for iPartition = CorrectionRange
                iPartition
                input_angles = AngleInfo(:,iPartition);
                input_kspace = kSpace(:,:,iPartition,:);
                ind_scan = find( input_angles == 0 );
                ScannedProj = 1:Nshot;
                if ~isempty(ind_scan) 
                    input_angles( ind_scan )        = [];
                    input_kspace(:,ind_scan,:,:)    = [];
                    ScannedProj( ind_scan )         = [];
                end

                Option.Indexes      = 1:length(input_angles);  
                Option.Angles       = input_angles;
                Option.InitAngle    = input_angles(1);
                Option.N            = 5;   % 1;
                Option.AngleRange   = '-G';
                if isfield( Option, 'PC_coeff' )
                    Option = rmfield( Option, 'PC_coeff' );
                end
                disp( 'Trajectory correction (RING)...' )
                for iCoil = 1:Ncoil
                    input_kspace2 = input_kspace(:,:,:,iCoil);
                    Option.PC_coeff = TrajectoryCorrection_RING( squeeze(input_kspace2), Option );
                    Option.PhaseCorrection = 'UsePC';
                    PCcorrected = RadialPhaseCorrection_RTCine( input_kspace2, input_kspace2, Option );
                    kSpace_TC(:,ScannedProj,iPartition,iCoil) = PCcorrected.kSpace;
                    clear PCcorrected
                end
            end
            clear  input_kspace  input_kspace2  input_angles
            Option.Angles = AngleInfo;
          


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% % FBP - density conpensation
fprintf( 'Calculate DCF for FBP...\n' )    
        Option.Nsamples         = Nsample;
        Option.AnglePrecision   = 4;
        Option.Display          = 0;
        Option.Angles           = AngleInfo;
        kSpaceEmpty_FBP         = abs(kSpace(:,:,:,1));
        kSpaceEmpty_FBP( find(kSpaceEmpty_FBP > 0) ) = 1;
        Option.WeightedContrast = kSpaceEmpty_FBP;
       
            DCF_FBP = gDCF_2DRadial( Option.Angles, Option );
            
     
        size(DCF_FBP)
% % %         figure(30); imagescn( abs(DCF_FBP),[0 1],[],[], 3); colormap jet


    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% % FBP - gridding
fprintf( 'Gridding radial k-space data...\n' )    
        ray_info_xy     = -0.5:1/Nsample:0.5-1/Nsample;
        ray_info_z      = -0.5:1/Npartition:0.5-1/Npartition;
        Grid_xy_FBP     = zeros( Nsample, Nshot, Npartition );
        Grid_z_FBP      = zeros( Nsample, Nshot, Npartition );
        for ii = 1:Npartition
            for jj = 1:Nshot
                Grid_xy_FBP(:,jj,ii) = ray_info_xy.*exp(1i*AngleInfo(jj,ii));
            end
            Grid_z_FBP(:,:,ii) = ray_info_z(ii);
        end
        Grid_xy_FBP     = 1i.*(-real(Grid_xy_FBP)) - imag(Grid_xy_FBP);
        Grid_z_FBP      = Grid_z_FBP*(-1);
        size(Grid_xy_FBP)
        size(Grid_z_FBP)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % manually exclude some bad coils
fprintf( 'Excluding bad coils...\n' )    
         if ~isempty( Index_ExcludingCoils )
             tmp_kspace = kSpace_TC;
             tmp_kspace(:,:,:,Index_ExcludingCoils) = [];
             kSpace_TC = tmp_kspace;
             clear tmp_kspace
         end
 [Nsample, Nshot, Npartition,  NumOfCoils] = size( kSpace_TC );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% % FBP using NUFFT
fprintf( 'NUFFT for FBP...\n' )    
        clear NUFFT_FBP
        NUFFT_FBP = MCNUFFT_3D_GPU_kxyzt_single( Grid_xy_FBP, Grid_z_FBP, DCF_FBP, ones(Nsample,Nsample,Npartition,NumOfCoils) );
        FBP_CoilSens = NUFFT_FBP'*kSpace_TC;
        
            FBP = sqrt( sum( abs(FBP_CoilSens).^2, 4 ) );
            FBP_norm = FBP./max(abs(FBP(:)));
            
            clear  FBP_norm  FBP   FBP_CoilSens_norm



        if strcmp( lower(DisplayFigure), 'on' )
            FBP = sqrt( sum( abs(FBP_CoilSens).^2, 4 ) );
            
            % Check this line for ZF recon
            FBP_norm = FBP./max(abs(FBP(:)));

            tmp_img = fliplr( imrotate( FBP_norm, RotatingAngle ) );
            [nx ny nz] = size( tmp_img );
            RangeX = nx/2+1-nx/4:nx/2+nx/4;
            MaxInt = 0.8
            figure(1); imagescn( abs(tmp_img), [0 MaxInt],[],[], 3 )
            figure(2); imagescn( abs(tmp_img(RangeX,RangeX,:)), [0 MaxInt],[],[], 3 )

            OtherView = permute( tmp_img(RangeX,RangeX,:), [1 3 2] );
            figure(3); imagescn( abs(OtherView), [0 MaxInt],[],[], 3 )
            OtherView = permute( tmp_img(RangeX,RangeX,:), [2 3 1] );
            figure(4); imagescn( abs(OtherView), [0 MaxInt],[],[], 3 )
            clear  FBP  FBP_norm  OtherView  FBP_CoilSens_norm
        end

        
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% % calculate 3D coil sens map
fprintf( 'Calculating coil sensitivity map...\n' )
        clear  CoilSens
        if 0
            % % for aniso-tropic resolution 
            parCoilSens.BlockSize = [ 6  6  4 ];
            [~, CoilSens] = adapt_array_2d_to_3d_pyo_BlockSize( FBP_CoilSens, parCoilSens );
        else
            % % for iso-tropic resolution 
            [~, CoilSens] = adapt_array_2d_to_3d_pyo( FBP_CoilSens );
            CoilSens = CoilSens./max(abs(CoilSens(:)));
        end
        
        clear tmp_nufft  dummy
        if strcmp( lower(DisplayFigure), 'on' )
            MaxCoilSignals = [];
            for iCoil = 1:size(kSpace_TC,4)
                tmp00 = kSpace_TC(:,:,:,iCoil);
                MaxCoilSignals(iCoil) = max(abs(tmp00(:)));
            end          
            MaxCoilSignals
            input_data = FBP_CoilSens;
            FBP_CoilSens_norm = input_data./max(abs(input_data(:)));
            clear  input_data  tmp00
            combined = [fliplr(imrotate(FBP_CoilSens_norm,RotatingAngle))   fliplr(imrotate(CoilSens,RotatingAngle))];
            sel_part = Npartition/2+1;
            figure(7); imagescn( abs(squeeze(combined(:,:,sel_part,:))),[0 0.99],[],[], 3 )
            sel_coil = 4   %NumOfCoils;
            figure(8); imagescn( abs(squeeze(combined(:,:,:,sel_coil))),[0 0.3],[],[], 3 )
            
% %             figure(9); imshow4( abs(combined)  )
        end
      
        clear NUFFT_FBP 
        clear FBP_CoilSens_norm   FBP_CoilSens  


 
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% % CS - rebinning k-space data based on respiratory motions
fprintf( 'Rebinning k-space data...\n' )
        if ~isempty( Index_ExcludingCoils_Navi )  && ( size(Navigator,4) == Ncoil )
            Navigator(:,:,:,Index_ExcludingCoils_Navi) = [];
            size( Navigator )
        end
        
        tmp_navi                                = Navigator(:,:,1:NaviSetting(2),:);
        param_SelfNavi.SelfNaviMethod           = NaviSetting(1);
        param_SelfNavi.TotalNumberOfNaviSpokes  = NaviSetting(2);  % # of Navigators
        param_SelfNavi.ActualNaviSpoke          =   NaviSetting(3);  %1=1st self-navi; 2=2nd; 3=3rd; 4=4th 
        param_SelfNavi.SpokeForDetectingBadNavi = NaviSetting(4);
        param_SelfNavi.CutMethod                = NaviSetting(5:6);
        param_SelfNavi.NumOfRespirationStates   = NaviSetting(7);
        Nrespiration                            = NaviSetting(7);
        param_SelfNavi.SelectPC                 =  SEL_PCA
        switch 2
            case 1,        param_SelfNavi.SignalAdjustPCA   = 'ON';    %'OFF';   %for TI chagnes along time
            case 2,        param_SelfNavi.SignalAdjustPCA   = 'OFF';   %for TI chagnes along time
        end
        if SEL_ImgNavi == 1
             param_SelfNavi.SelfNaviMethod = 7
        end
        switch SEL_VS
            case 0,     ViewSharing = 'NoVS';   ViewSharingParam = [ 0  0   0 ];
            case 1,     ViewSharing = 'WithVS'; ViewSharingParam = [ 0.3  3  round(Nsample*0.1) ];  %[ amount of VS from neighbors(0-1), # of excluding kz,  VS width ]
        end
        
        switch   param_SelfNavi.SelfNaviMethod
            case 6
                param_SelfNavi.Display      = DisplayFigure;
                param_SelfNavi.NaviWindow   = [1 Nsample];
                [Resp_Signal, ExtraInfo]    = SelfNavigation_StackOfStars_v6( tmp_navi,  param_SelfNavi );

                Index_DataSelection = ExtraInfo.Index_KeepData;
                NumOfIntervals      = ExtraInfo.NumOfIntervals;

            case 7
            % % for 2D image navi...
            % %  : NOT USED
                SearchRange                             = 201:300;
                param_SelfNavi.Display                  = DisplayFigure;
                param_SelfNavi.ExcludingBadSignal       =  'ON';  %'ON'; 'OFF'
                param_SelfNavi.NaviWindow               = [1 Nsample];
                param_SelfNavi.SearchRange              = SearchRange;
                param_SelfNavi.Nrespiration             = Nrespiration;
                param_SelfNavi.NaviSignal_1D_OnLiver    = NaviSignal_1D_OnLiver;
                switch  2
                    case 1,     SelectRebinMethod = 'PCA'
                    case 2,     SelectRebinMethod = 'ImgNaviPCA'
                    case 3,     SelectRebinMethod = 'ImgNaviDiff'
                    case 4,     SelectRebinMethod = 'ImgNaviEFA'
                    otherwise,  SelectRebinMethod = 'PCA'
                end
                param_SelfNavi.SelectRebinMethod    = SelectRebinMethod;
                Result_SelfNavi                     = SelfNavigation_StackOfStars_v7( tmp_navi,  param_SelfNavi );

                Index_DataSelection                 = Result_SelfNavi.Index_KeepData;
                Resp_Signal                         = Result_SelfNavi.Resp_Signal_ImgNavi;
                
                SortedNaviValues                    = Result_SelfNavi.SortedNaviValues;
                RespRebin_RefIndex                  = Result_SelfNavi.SortedNaviIndices;
                param_Rebinning.InitMaxRespBin      = 2;
                param_Rebinning.RebinWindow         = Result_SelfNavi.IncStep; 
                param_Rebinning.ViewSharing         = ViewSharing; 
                param_Rebinning.ViewSharingParam    = ViewSharingParam; 
                param_Rebinning.Display             = DisplayFigure;
                Rebinned_Info                       = AutoRebinning_Respiration( SortedNaviValues, RespRebin_RefIndex, param_Rebinning );
            otherwise
                % do nothing
        end
        if  SEL_ImgNavi == 0  
            if param_SelfNavi.SelfNaviMethod == 6
                NaviInfo = sprintf( 'M%dNavi%dEFA%d', param_SelfNavi.SelfNaviMethod, param_SelfNavi.ActualNaviSpoke, 2 )  %param_SelfNavi.SelectPC )
            else
                NaviInfo = sprintf( 'M%dNavi%dPCA%d', param_SelfNavi.SelfNaviMethod, param_SelfNavi.ActualNaviSpoke, param_SelfNavi.SelectPC )
            end
        else
            NaviInfo = sprintf( 'M%d%s', param_SelfNavi.SelfNaviMethod, SelectRebinMethod )
            if SEL_VS == 1
                if strcmp( SelectRebinMethod, 'ImgNaviPCA' ) 
                   VS_amount = sprintf( 'VS%1.1f',  ViewSharingParam(1) );
                else
                   VS_amount = sprintf( 'VS' );
                end
            end
        end
        
        
    if  SEL_ImgNavi == 0   %isempty( Rebinned_Info )
        [B, RespRebin_RefIndex] = sort( Resp_Signal, 'ascend' );
        MaxNavi = max( B(:) );
        MinNavi = min( B(:) );
        IncStep = min( 1/Nrespiration,  1/NumOfIntervals );
        B_norm = (B - MinNavi)./(MaxNavi - MinNavi);
        Ranges = 0:IncStep:1;
        if strcmp( lower(DisplayFigure), 'on' )
            figure%(555); 
            subplot(1,2,1); plot( Resp_Signal ); axis square;
            subplot(1,2,2); plot( B_norm ); 
            for ii = 1:length(Ranges), hold on; plot( 1:length(B), Ranges(ii).*ones(1,length(B)) ); end;  axis square; 
        end
        
        IndexRespBins   = {};
        val_navi        = {};
        MaxSpokes       = 0;
        MaxRespBin      = 0;
        ExcludingCells  = [];
        CombineInd      = [];
        for ii = 1:length(Ranges)-1
            [tmp_ind,tmp_val ]  = find( B_norm >= Ranges(ii) & B_norm <= Ranges(ii+1) );
            IndexRespBins{ii}   = RespRebin_RefIndex( tmp_ind );
            val_navi{ii}        = B_norm( tmp_ind );
            if MaxSpokes < length(IndexRespBins{ii}), MaxSpokes = length(IndexRespBins{ii}); MaxRespBin = ii; end
        end
        IndexRespBins
        Init_MaxRespBin = MaxRespBin;

        param_Rebinning.InitMaxRespBin  = Init_MaxRespBin;
        param_Rebinning.RebinWindow     = IncStep; 
        param_Rebinning.Display         = 'on'  %DisplayFigure;
        [Rebinned_Info]                 = AutoRebinning_Respiration( B_norm, RespRebin_RefIndex, param_Rebinning );
    end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% % k-space rebinning based on the Rebinned_Info
    Rebinned_Info
            % Number of bins:  [ below,   MaxResp,  above];
            Index_NumOfRespBins = [ 1 1 1 ];
            tmp_length_below = length(Rebinned_Info.Below)
            if tmp_length_below < 60 
                Index_NumOfRespBins(1) = 0;
            else
                if tmp_length_below >= 150  && tmp_length_below < 300
                    Index_NumOfRespBins(1) = 2;
                elseif tmp_length_below >= 300  && tmp_length_below < 400
                    Index_NumOfRespBins(1) = 3;
                elseif tmp_length_below >= 400  && tmp_length_below < 500
                    Index_NumOfRespBins(1) = 4;
                elseif tmp_length_below >= 500  %&& tmp_length_below < 550
                    Index_NumOfRespBins(1) = 5;
                end
            end
            tmp_length_above = length(Rebinned_Info.Above)
            if tmp_length_above < 60 
                Index_NumOfRespBins(3) = 0;
            else
                if tmp_length_above >= 150  && tmp_length_above < 300
                    Index_NumOfRespBins(3) = 2;
                elseif tmp_length_above >= 300  && tmp_length_above < 400
                    Index_NumOfRespBins(3) = 3;
                elseif tmp_length_above >= 400  && tmp_length_above < 500
                    Index_NumOfRespBins(3) = 4;
                elseif tmp_length_above >= 500  %&& tmp_length_above < 550
                    Index_NumOfRespBins(3) = 5;
                end
            end
            Index_NumOfRespBins


            IndexRespBins   = {};
            IndBelow        = Index_NumOfRespBins(1);
            IndMaxResp      = Index_NumOfRespBins(2);
            IndAbove        = Index_NumOfRespBins(3);
            if IndBelow == 1
                IndexRespBins{1} = RespRebin_RefIndex( Rebinned_Info.Below );
            elseif IndBelow > 1
                tmp_length = length( Rebinned_Info.Below );
                tmp_jump = 1:round(tmp_length/IndBelow):tmp_length;
                for iBin = 1:IndBelow-1
                    IndexRespBins{iBin} = RespRebin_RefIndex( Rebinned_Info.Below(1:tmp_jump(iBin+1)) );
                end
                IndexRespBins{IndBelow} = RespRebin_RefIndex( Rebinned_Info.Below(tmp_jump(iBin+1)+1:end) );
            end
            IndAbove_begin = max( IndBelow+IndMaxResp+1, 3 );
            if IndAbove == 1
                IndexRespBins{IndAbove_begin} = RespRebin_RefIndex( Rebinned_Info.Above );
            elseif IndAbove > 1
                tmp_length = length( Rebinned_Info.Above );
                tmp_jump = 1:round(tmp_length/IndAbove):tmp_length;
                for iBin = 1:IndAbove-1
                    IndexRespBins{IndAbove_begin+iBin-1} = RespRebin_RefIndex( Rebinned_Info.Above(tmp_jump(iBin):tmp_jump(iBin+1)) );
                end
                IndexRespBins{IndAbove_begin+IndAbove-1} = RespRebin_RefIndex( Rebinned_Info.Above(tmp_jump(iBin+1)+1:end) );
            end
            if IndBelow == 0
                IndexRespBins{1} = IndexRespBins{end};
                IndexRespBins(end) = [];
            end
            if SEL_VS == 1
                MaxRespIndex = Rebinned_Info.MaxRespIndex_VS;
            else
                MaxRespIndex = Rebinned_Info.MaxRespIndex;
            end
            tm_MaxRespBinNumber = max( IndBelow+IndMaxResp, 2 );
            IndexRespBins{tm_MaxRespBinNumber} = RespRebin_RefIndex( MaxRespIndex );
            MaxRespBin = tm_MaxRespBinNumber;

            tmp_Nrespiration = size(IndexRespBins,2);
            MaxSpokes = size(  IndexRespBins{MaxRespBin}, 1 );
            IndexRespBins
            [ MaxSpokes  MaxRespBin ]


        kSpace_adjusted     = kSpace_TC( :, Index_DataSelection, :, : );
        AngleInfo_adjusted  = AngleInfo( Index_DataSelection, : );
        kSpace_Rebin        = zeros( Nsample, MaxSpokes, Npartition, NumOfCoils, tmp_Nrespiration );
        Angle_Rebin         = zeros( MaxSpokes, Npartition, tmp_Nrespiration );
        for iResp = 1:tmp_Nrespiration
            ind = IndexRespBins{iResp};
            kSpace_Rebin(:,1:length(ind),:,:,iResp) = kSpace_adjusted(:,ind,:,:);
            Angle_Rebin(1:length(ind),:,iResp) = AngleInfo_adjusted(ind,:);
        end
        kSpaceEmpty_CS = abs(squeeze(kSpace_Rebin(:,:,:,1,:)));
        kSpaceEmpty_CS( find( kSpaceEmpty_CS > 0 ) ) = 1;
        [NX NY NZ NC Nresp] = size(kSpace_Rebin);
        clear kSpace_adjusted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    



    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% % CS - density conpensation
fprintf( 'Calculating DCF for CS...\n' )
        Option.Nsamples         = Nsample;
        Option.AnglePrecision   = 4;
        Option.Display          = 0;
        DCF_CS = [];
        for iResp = 1:Nresp
            iResp
            if SEL_VS == 1
                if iResp == MaxRespBin
                    VS_info = [ length(Rebinned_Info.MaxRespIndex_VS_before),   length(Rebinned_Info.MaxRespIndex_VS_after) ];
                    VS_kz = ViewSharingParam(2);
                    VS_Width = ViewSharingParam(3);
                    Mask_VS = ones( size(kSpaceEmpty_CS(:,:,:,iResp)) );
                    for iKz = Npartition/2+1-VS_kz:Npartition/2+1+VS_kz
                        if VS_info(1) > 0
                            Mask_VS( Nsample/2+1-VS_Width:Nsample/2+1+VS_Width, 1:VS_info(1), iKz ) = 0;
                        end
                        Mask_VS( Nsample/2+1-VS_Width:Nsample/2+1+VS_Width, end+1-VS_info(2):end, iKz ) = 0;
                    end
                    kSpaceEmpty_CS(:,:,:,iResp) = kSpaceEmpty_CS(:,:,:,iResp).*Mask_VS;
                end
                VS_amount = sprintf( 'VS%1.1f',  ViewSharingParam(1) );
                 VS_kz = ViewSharingParam(2);
                VS_Width = ViewSharingParam(3);
                VS_shape = sprintf( '%s%02d%03d', VS_amount, VS_kz, VS_Width );
               
            else
                VS_shape = 'none';
            end
            Option.WeightedContrast = kSpaceEmpty_CS(:,:,:,iResp);
            tmp_dcf = gDCF_2DRadial( Angle_Rebin(:,:,iResp), Option );
            DCF_CS(:,:,:,iResp) = tmp_dcf; 
        end
        clear tmp_dcf
% %     size(DCF_CS)
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% % CS - gridding
fprintf( 'Gridding CS k-space data...\n' )
        ray_info_xy     = -0.5:1/Nsample:0.5-1/Nsample;
        ray_info_z      = -0.5:1/NZ:0.5-1/NZ;
        Grid_xy_CS      = zeros( Nsample, NY, NZ, Nresp );
        Grid_z_CS       = zeros( Nsample, NY, NZ );
        for ii = 1:Nresp
            for jj = 1:NZ
                for kk = 1:NY
                    Grid_xy_CS(:,kk,jj,ii) = ray_info_xy.*exp( 1i*Angle_Rebin(kk,jj,ii) );
                end
                if ii == 1
                    Grid_z_CS(:,:,jj) = ray_info_z(jj);
                end
            end
        end
        Grid_xy_CS  = 1i.*(-real(Grid_xy_CS)) - imag(Grid_xy_CS);
        Grid_z_CS   = Grid_z_CS*(-1);
% % %     size(Grid_xy_CS)
% % %     size(Grid_z_CS)
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% % CS - NUFFT and preparing parameters
fprintf( 'Preparing CS reconstruction...\n' )
        clear param
        param.E = MCNUFFT_3D_GPU_kxyzt_kSpaceEmpty( Grid_xy_CS, Grid_z_CS, DCF_CS, CoilSens, kSpaceEmpty_CS );
        param.y = kSpace_Rebin;
        Recon_ZeroFilled = param.E'*param.y;
        size(Recon_ZeroFilled)
        if strcmp( lower(DisplayFigure), 'on' )
            img_norm = fliplr(imrotate(Recon_ZeroFilled,RotatingAngle))./max(abs(Recon_ZeroFilled(:)));
            [NX3 NY3 NZ3 Nresp] = size(img_norm);
            sel_resp = MaxRespBin
            OverSamplingFactor = NZ3/1.2; % 20% oversampling in kz
            sel_z = NZ3/2+1-OverSamplingFactor/2:NZ3/2+OverSamplingFactor/2;
            MaxInt = 0.99
            figure(11); imagescn( abs(img_norm(:,:,sel_z,sel_resp)),[0 MaxInt],[],[], 3 )
            figure(12); imagescn( abs(squeeze(img_norm(:,:,Npartition/2+1,:))),[0 MaxInt],[],[], 3 )
    % %         figure(13); imagescn( abs(squeeze(img_norm(:,:,10,:))),[0 MaxInt],[],[], 3 )
            clear img_norm
        end
    
        MaxProjNumber = size( IndexRespBins{MaxRespBin}, 1  )

        if SEL_SaveData == 1
            tmp_name    = FileList{Sel_File}(1:end-4);
            tmp_ind     = strfind( tmp_name, 'FID' );
            subname     = sprintf( 'Recon_%s', tmp_name(tmp_ind:end) );
            if SEL_VS == 1
                filename2 = [path_SaveData subname sprintf( '_NX%03d_Proj%03d_Part%03d_Resp%03d_Coil%02d_%s_%s__Iter00.mat', ...
                    NX/2, MaxProjNumber, NZ, param_SelfNavi.NumOfRespirationStates, NC, NaviInfo, VS_shape ) ]
            else
                filename2 = [path_SaveData subname sprintf( '_NX%03d_Proj%03d_Part%03d_Resp%03d_Coil%02d_%s__Iter00.mat', ...
                    NX/2, MaxProjNumber, NZ, param_SelfNavi.NumOfRespirationStates, NC, NaviInfo ) ]
            end
            save( filename2, 'Recon_ZeroFilled', 'MaxRespBin', 'IndexRespBins', '-v7.3' );
        end

    
        param.TPCA              = TempPCA_3D();
        param.TTV               = TV_3D_Temp();
        param.STV               = TV_3DPERF();
        param.nite              = 8;
        param.display           = 1;
        param.TPCAw             = 0;
        param.STVw              = 0;
        param.VW                = 'on';
        param.gpuNumber         = sel_gpu;
        param.NaviInfo          =  NaviInfo;
        param.ViewSharingInfo   =  VS_shape
        ITER                    = 3;

    
        clear kSpace_Rebin    kSpaceEmpty_CS  
% % %         clear  CoilSens   kSpace   kSpace_TC
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% % CS - iterative recon
fprintf( 'Iterative CS reconstruction...\n' )
    
        for lambda =  0.001   
            param.Lambda    = lambda;
            param.TTVw      = lambda;
            param.TPCAw     = 0;  % param.TTVw/2;
            param

            Recon_CS_TTV_VW = Recon_ZeroFilled;
            for ite = 1:ITER
                fprintf( 'Iteration loop: %d out of %d\n', ite, ITER )
                tic
                Recon_CS_TTV_VW = CSL1NlCg_TTV_TPCA_STV_VW( Recon_CS_TTV_VW, param );
                toc
            end
            if SEL_SaveData == 1
                if SEL_VS == 1
                    filename_save = [path_SaveData  subname  sprintf( '_NX%03d_Proj%03d_Part%03d_Resp%03d_Coil%02d_TTV%s_TPCA%s_%s_%s_Iter%02d.mat', ...
                        NX/2, MaxProjNumber, NZ, param_SelfNavi.NumOfRespirationStates, NC, sprintf( '%1.6f', param.TTVw ), num2str(param.TPCAw), ...
                        NaviInfo, VS_shape, ite*(param.nite+1) )  ]
                else
                    filename_save = [path_SaveData  subname  sprintf( '_NX%03d_Proj%03d_Part%03d_Resp%03d_Coil%02d_TTV%s_TPCA%s_VW%s_%s_Iter%02d.mat', ...
                        NX/2, MaxProjNumber, NZ, param_SelfNavi.NumOfRespirationStates, NC, sprintf( '%1.6f', param.TTVw ), num2str(param.TPCAw), param.VW, NaviInfo, ite*(param.nite+1) )  ]
                end
                save( filename_save, 'Recon_CS_TTV_VW', 'MaxRespBin', 'IndexRespBins', '-v7.3' );
            end
        end % lambda

            CS_recon            = fliplr(imrotate(Recon_CS_TTV_VW,RotatingAngle))./max(abs(Recon_CS_TTV_VW(:)));
            [nx ny nz nt]       = size( CS_recon );
            RangeX              = nx/2+1-nx/4:nx/2+nx/4;
            OverSamplingFactor  = nz/1.2; % 20% oversampling in kz
            sel_z               = nz/2+1-OverSamplingFactor/2:nz/2+OverSamplingFactor/2;
            Img_CSrecon         = CS_recon(RangeX,RangeX,sel_z,:);
            
            ParamLRBW.block_size    = [ 6  6  6 ];  %[4 4 4];
            ParamLRBW.tau_size      = 0.5;
            ParamLRBW.Iterations    = 1;
            Img_CSrecon_lrbw        = LRBW_Denoising_3Dt( Img_CSrecon, ParamLRBW ); 
            
            MaxInt      = 0.99;
            combined    = [ Img_CSrecon   Img_CSrecon_lrbw ];
            figure(201); imagescn( abs(combined(:,:,:,MaxRespBin)), [0  MaxInt],[],[],3 )
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % the final result of CS reconstruction:
            
            img = abs(Img_CSrecon_lrbw(:,:,:,MaxRespBin));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

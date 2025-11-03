function Report = SecondaryCapture_to_NII(examdir,all_masks_on_same_series)
    % When segmentations of an image are saved as secondary captures, the 
    % resultant DICOM files contain no position or orientation information.
    % Further, it is often the case that files possess no information about
    % which image series was used as the basis of the segmentation. The
    % secondary captures are often saved after some degree of zooming
    % and/or panning have been applied. Lastly, it is often the case that
    % secondary captures are only saved for a subset of slices of the
    % orginal image.
    % 
    % So here is a function that finds which original image series was used
    % to draw a given segmentation and creates a nifti file of the
    % segmentation
    %
    % inputs:
    %   examdir - full file path of the exam directory containing the DICOM
    %   images and secondary captures (e.g., */PT-*/ST* or
    %   */Patient^Name/yyyyMMdd)

    arguments
        examdir {mustBeFolder}
        all_masks_on_same_series {mustBeMember(all_masks_on_same_series,[0 1])} = 1
    end

    % First, find all secondary captures
    dcms = dir(fullfile(examdir,'*','*.dcm'));
    % ignore structured reports, presentation modes (these files, I think,
    % are always named SR*.dcm and PR*.dcm, respectively
    dcms(cellfun(@(f) ismember(f(1:2),{'SR','PR'}),{dcms.name})) = [];

    % Now, read the ImageType field from all dicoms in each folder (the
    % reason we read all dicoms instead of just the first is because I have
    % seen cases where primary and secondary captures are mixed in the same
    % scan folder.
    scans = unique({dcms.folder});
    metadata = cellfun(@(f) fast_dcm_head(f,{'ImageType'}),scans,'Uni',0);
    is_secondary = cellfun(@(f) all(contains({f.ImageType},'secondary','IgnoreCase',true)),metadata);
    is_primary = cellfun(@(f) all(contains({f.ImageType},'primary','IgnoreCase',true)),metadata);
    is_volumetric = cellfun(@length,metadata)>1;
    secondary_scans = scans(is_volumetric & is_secondary);
    primary_scans = scans(is_volumetric & is_primary);

    % Now, convert each primary scan to nifti, read
    niidir = fullfile(examdir,'NII');
    cellfun(@(f) dcm2niix(f,niidir),primary_scans,'Uni',0);
    list = dir(fullfile(niidir,'*.nii.gz'));
    list = fullfile({list.folder},{list.name});
    orig = struct();
    not3d = false(size(list));
    for n=1:length(list)
        orig(n).list = list{n};
        orig(n).head = niftiinfo(list{n});
        if numel(orig(n).head.ImageSize)~=3
            not3d(n) = true;
        elseif orig(n).head.ImageSize(3)==1
            not3d(n) = true;
        end
    end
    % Presumably not looking at 4D images (but who knows, right?!)
    orig(not3d) = [];

    if all(strncmp({dcms.name},'CT',2))
        % Modality is CT, read in images, apply scaling + offset, apply
        % windowing
        parfor n=1:length(orig)
            temp = double(niftiread(orig(n).list));
            temp = temp.*orig(n).head.MultiplicativeScaling + orig(n).head.AdditiveOffset;
            orig(n).Primary_I = mat2gray(temp,[-160 240]);
        end
    else
        parfor n=1:length(orig)
            orig(n).Primary_I = imadjustn(mat2gray(niftiread(orig(n).list)));
        end
    end

    if all_masks_on_same_series
        fprintf('Looking for matches for %s\n',strrep(examdir,'\','\\'));
        [Secondary_I,~,~,Secondary_Masks] = combine_secondaryCaptures(secondary_scans);
        % instantiate output
        Report = table({Secondary_Masks.list}','VariableNames',{'Secondary_Image'});
        % Get "Speeded Up Robust Features" from the image
        parfor slice=1:size(Secondary_I,3)
            surfpoints = detectSURFFeatures(Secondary_I(:,:,slice),'NumScaleLevels',6);
            f1{slice} = extractFeatures(Secondary_I(:,:,slice),surfpoints);
            % Also get them flipped in case the capture was drawn on a
            % flipped image
            surfpoints = detectSURFFeatures(flip(Secondary_I(:,:,slice)),'NumScaleLevels',6);
            f1_flip{slice} = extractFeatures(flip(Secondary_I(:,:,slice)),surfpoints);
        end
        % Create a copy for later
        si_copy = Secondary_I;

        % Quickly sort potential matches 
        Iz = cellfun(@(f) squeeze(mean(mean(f))),{orig.Primary_I},'Uni',0);
        Sz = squeeze(mean(mean(Secondary_I,'omitnan'),'omitnan'));
        maxcorr = zeros(length(Iz),1);
        for p=1:length(Iz)
            maxcorr(p) = quick_match(Sz,Iz{p});
        end
        [~,order] = sort(maxcorr,'descend');
        orig = orig(order);
        
        rho = NaN(length(orig),1);
        flipdim1 = false(length(orig),1);
        for k=1:length(orig)
            % Reset secondary capture
            Secondary_I = si_copy;
            MIP_Secondary = max(Secondary_I,[],3);
    
            % Read in the image from origlist
            [~,on] = fileparts(orig(k).list);
            fprintf(sprintf('\t%s: ',on));
            I = orig(k).Primary_I;
    
            % Extract features
            f2 = cell(0);
            parfor slice=1:size(I,3)
                surfpoints = detectSURFFeatures(I(:,:,slice),'NumScaleLevels',6);
                f2{slice} = extractFeatures(I(:,:,slice),surfpoints);
            end
            
            % Find the substack that optimizes the number of matched points
            [cI,orig(k).idx,corresponds] = find_substack(f2,f1,I);
            if ~corresponds
                % No correspondance, try flipping
                flipdim1(k) = true;
                Secondary_I = flip(Secondary_I,1);
                MIP_Secondary = max(Secondary_I,[],3);
                [cI,orig(k).idx,corresponds] = find_substack(f2,f1_flip,I);
                if ~corresponds
                    fprintf('No slicewise correspondance.\n')
                    continue
                end
            end
            orig(k).refview = imref3d(size(cI));
    
            % Create maximum intensity projection of original image series
            MIP_orig = max(cI,[],3);
    
            tform = matchedpoint_registration(MIP_Secondary,MIP_orig);
            try
                rMIP_Secondary = imwarp(MIP_Secondary,tform,'OutputView',imref2d(size(MIP_orig)));
                rMIP_Secondary(rMIP_Secondary==0) = NaN;
            catch
                rMIP_Secondary = NaN(size(MIP_orig));
            end
            if isempty(tform) || corr(MIP_orig(:),rMIP_Secondary(:),'rows','complete')<0.75
                flipdim1(k) = ~flipdim1(k);
                Secondary_I = flip(Secondary_I,1);
                MIP_Secondary = max(Secondary_I,[],3);
                [cI,orig(k).idx,corresponds] = find_substack(f2,f1_flip,I);
                if ~corresponds
                    fprintf('No slicewise correspondance.\n')
                    continue
                end
                orig(k).refview = imref3d(size(cI));
                % Create maximum intensity projection of original image series
                MIP_orig = max(cI,[],3);
                tform = matchedpoint_registration(MIP_Secondary,MIP_orig);
                try
                    rMIP_Secondary = imwarp(MIP_Secondary,tform,'OutputView',imref2d(size(MIP_orig)));
                    rMIP_Secondary(rMIP_Secondary==0) = NaN;
                catch
                    rMIP_Secondary = NaN(size(MIP_orig));
                end
                if isempty(tform) || corr(MIP_orig(:),rMIP_Secondary(:),'rows','complete')<0.75
                    fprintf('could not estimate transform.\n')
                    continue
                end
            end
    
            % Round to nearest 90 degree rotation and apply warp
            tform.R = round(tform.R); % This adjusts the rotation and automatically applies to the similarity matrix A
            mat = tform.A;
            mat = [mat(1:2,:);0 0 0;mat(3,:)];
            mat = [mat(:,1:2) [0;0;1;0] mat(:,3)];
            orig(k).tform3d = affinetform3d(mat);
            Secondary_I(round(0.9*size(Secondary_I,1)):end,round(0.75*size(Secondary_I,2)):end,:) = NaN;
            temp = imwarp(Secondary_I,orig(k).tform3d,'OutputView',orig(k).refview);
            temp(temp==0) = NaN;
            cI(cI==0) = NaN;
            rho(k) = corr(temp(:),cI(:),'rows','complete');
            fprintf('r = %0.3f\n',rho(k))
            if rho(k)>0.97
                % This is either the right match, or close enough
                break
            end
        end

        % Now use the best match according to correlation coefficients
        [maxrho,irho] = max(rho);
        if ~isnan(maxrho)
            for m=1:length(Secondary_Masks)
                Report.Max_Image_Correlation(m) = maxrho;
                Report.IsMatched(m) = true;
                % Transform Secondary Capture to Original Image
                if flipdim1(irho)
                    Secondary_Masks(m).BW = flip(Secondary_Masks(m).BW);
                end
                mask_nii = write_matched_mask(orig(irho),Secondary_Masks(m));
                Report.Matched_Image{m} = orig(irho).list;
                Report.Mask_Name{m} = mask_nii;
            end
            
        end

    else
        % instantiate output
        Report = table(secondary_scans','VariableNames',{'Secondary_Image'});
        % For each of the secondary scans, find the matching image
        for sc=1:length(secondary_scans)
            [~,secnam] = fileparts(secondary_scans{sc});
            fprintf('Looking for matches for %s\n',secnam);
            % Get the dicom file metadata
            mask_dcms = struct2table(fast_dcm_head(secondary_scans{sc},{'InstanceNumber','Rows','Columns'}));
            mask_dcms = sortrows(mask_dcms,'InstanceNumber');
            % Some directories are "dirty" with combined secondary and original
            % images. Toss out any with lower pixel count. Why lower? It worked
            % for the use-case this was developed for...it might not be the
            % right call
            mask_dcms(mask_dcms.Rows<max(mask_dcms.Rows),:) = [];
            % Keep only unique instance numbers (for rare cases when secondary
            % captures have apparently been re-saved in the same directory).
            [~,ia] = unique(mask_dcms.InstanceNumber);
            mask_dcms = mask_dcms(ia,:);
            
            if mask_dcms.Rows(1)~=mask_dcms.Columns(1)
                % non-square images (which are often volumetric renderings)
                continue
            end
            % Instantiate arrary
            M = uint8(zeros([mask_dcms.Rows(1) mask_dcms.Columns(1) 3 height(mask_dcms)]));
            for slice=1:height(mask_dcms)
                % Get the pixel data
                M(:,:,:,slice) = dicomread(mask_dcms.Filename{slice});
            end
            % swap 3rd and 4th dimensions to put color channels in dim4
            M = permute(M,[1 2 4 3]);
        
            % compute correlation between channels
            chancorr = corr(double(reshape(M,[],3)));
            % get the index of the channel unlike the other two
            diffi = find(sum(chancorr>0.999)==1);
            % ensure mask channel is 1 (red)
            if diffi==2
                M = M(:,:,:,[2 1 3]);
            elseif diffi==3
                M = M(:,:,:,[3 2 1]);
            end
    
            % Get the grayscale image
            Secondary_I = mat2gray(M(:,:,:,2));
            % replace the mask area with NaN
            Secondary_I(M(:,:,:,1)~=M(:,:,:,2)) = NaN;
            % Get the mask
            Secondary_Mask = (double(M(:,:,:,1))./double(M(:,:,:,2)))>=1.5;
            % Get rid of very small, spurious clusters in the image
            Secondary_Mask = ~bwareaopen(~Secondary_Mask,5,6);
    
            % Get "Speeded Up Robust Features" from the image
            parfor slice=1:size(Secondary_I,3)
                surfpoints = detectSURFFeatures(Secondary_I(:,:,slice),'NumScaleLevels',6);
                f1{slice} = extractFeatures(Secondary_I(:,:,slice),surfpoints);
                % Also get them flipped in case the capture was drawn on a
                % flipped image
                surfpoints = detectSURFFeatures(flip(Secondary_I(:,:,slice)),'NumScaleLevels',6);
                f1_flip{slice} = extractFeatures(flip(Secondary_I(:,:,slice)),surfpoints);
            end
            % Create a copy for later
            si_copy = Secondary_I;
    
            % Quickly sort potential matches 
            Iz = cellfun(@(f) squeeze(mean(mean(f))),Primary_I,'Uni',0);
            Sz = squeeze(mean(mean(Secondary_I,'omitnan'),'omitnan'));
            maxcorr = zeros(length(Iz),1);
            for p=1:length(Iz)
                maxcorr(p) = quick_match(Sz,Iz{p});
            end
            [~,order] = sort(maxcorr,'descend');
            origlist = origlist(order);
            orighead = orighead(order);
            Primary_I = Primary_I(order);
    
            
            rho = NaN(length(origlist),1);
            tform3d = cell(length(origlist),1);
            flipdim1 = false(length(origlist),1);
            refview = cell(length(origlist),1);
            idx = cell(length(origlist),1);
            for k=1:length(origlist)
                % Reset secondary capture
                Secondary_I = si_copy;
                MIP_Secondary = max(Secondary_I,[],3);
        
                % Read in the image from origlist
                [~,on] = fileparts(origlist{k});
                fprintf(sprintf('\t%s: ',on));
                I = Primary_I{k};
        
                % Extract features
                f2 = cell(0);
                parfor slice=1:size(I,3)
                    surfpoints = detectSURFFeatures(I(:,:,slice),'NumScaleLevels',6);
                    f2{slice} = extractFeatures(I(:,:,slice),surfpoints);
                end
                
                % Find the substack that optimizes the number of matched points
                [cI,idx{k},corresponds] = find_substack(f2,f1,I);
                if ~corresponds
                    % No correspondance, try flipping
                    flipdim1(k) = true;
                    Secondary_I = flip(Secondary_I,1);
                    MIP_Secondary = max(Secondary_I,[],3);
                    [cI,idx{k},corresponds] = find_substack(f2,f1_flip,I);
                    if ~corresponds
                        fprintf('No slicewise correspondance.\n')
                        continue
                    end
                end
                refview{k} = imref3d(size(cI));
        
                % Create maximum intensity projection of original image series
                MIP_orig = max(cI,[],3);
        
                tform = matchedpoint_registration(MIP_Secondary,MIP_orig);
                try
                    rMIP_Secondary = imwarp(MIP_Secondary,tform,'OutputView',imref2d(size(MIP_orig)));
                    rMIP_Secondary(rMIP_Secondary==0) = NaN;
                catch
                    rMIP_Secondary = NaN(size(MIP_orig));
                end
                if isempty(tform) || corr(MIP_orig(:),rMIP_Secondary(:),'rows','complete')<0.75
                    flipdim1(k) = ~flipdim1(k);
                    Secondary_I = flip(Secondary_I,1);
                    MIP_Secondary = max(Secondary_I,[],3);
                    [cI,idx{k},corresponds] = find_substack(f2,f1_flip,I);
                    if ~corresponds
                        fprintf('No slicewise correspondance.\n')
                        continue
                    end
                    refview{k} = imref3d(size(cI));
                    % Create maximum intensity projection of original image series
                    MIP_orig = max(cI,[],3);
                    tform = matchedpoint_registration(MIP_Secondary,MIP_orig);
                    try
                        rMIP_Secondary = imwarp(MIP_Secondary,tform,'OutputView',imref2d(size(MIP_orig)));
                        rMIP_Secondary(rMIP_Secondary==0) = NaN;
                    catch
                        rMIP_Secondary = NaN(size(MIP_orig));
                    end
                    if isempty(tform) || corr(MIP_orig(:),rMIP_Secondary(:),'rows','complete')<0.75
                        fprintf('could not estimate transform.\n')
                        continue
                    end
                end
        
                % Round to nearest 90 degree rotation and apply warp
                tform.R = round(tform.R); % This adjusts the rotation and automatically applies to the similarity matrix A
                mat = tform.A;
                mat = [mat(1:2,:);0 0 0;mat(3,:)];
                mat = [mat(:,1:2) [0;0;1;0] mat(:,3)];
                tform3d{k} = affinetform3d(mat);
                Secondary_I(round(0.9*size(Secondary_I,1)):end,round(0.75*size(Secondary_I,2)):end,:) = NaN;
                temp = imwarp(Secondary_I,tform3d{k},'OutputView',refview{k});
                temp(temp==0) = NaN;
                cI(cI==0) = NaN;
                rho(k) = corr(temp(:),cI(:),'rows','complete');
                fprintf('r = %0.3f\n',rho(k))
                if rho(k)>0.97
                    % This is either the right match, or close enough
                    break
                end
            end
    
            % Now use the best match according to correlation coefficients
            [maxrho,irho] = max(rho);
            Report.Max_Image_Correlation(sc) = maxrho;
            if ~isnan(maxrho)
                Report.IsMatched(sc) = true;
                % Transform Secondary Capture to Original Image
                if flipdim1(irho)
                    Secondary_Mask = flip(Secondary_Mask);
                end
                
                Report.Matched_Image{sc} = origlist{irho};
                Report.Mask_Name{sc} = [mname '.gz'];
            else
                Report.IsMatched(sc) = false;
            end
    
    
        end
    end


    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    % Subroutines
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    function mname = write_matched_mask(O,S)
        alignedMask = uint8(zeros(O.head.ImageSize));
        tempmsk = imwarp(S.BW,O.tform3d,'OutputView',O.refview);
        alignedMask(:,:,O.idx) = tempmsk;
        O.head.Datatype = 'uint8';
        O.head.MultiplicativeScaling = 1;
        O.head.AdditiveOffset = 0;
        [~,se] = fileparts(S.list);
        mname = strrep(O.list,'.nii.gz',['_' se '.nii']);
        niftiwrite(alignedMask,mname,O.head)
        switch O.head.TransformName
            case 'Sform'
                % Add the qform
                fid = fopen(mname,'r+');
                fseek(fid,76,'bof');
                fwrite(fid,O.head.raw.pixdim(1),'float');
                fseek(fid,252,'bof');
                fwrite(fid,O.head.raw.qform_code,'short');
                fseek(fid,256,'bof');
                fwrite(fid,O.head.raw.quatern_b,'float');
                fwrite(fid,O.head.raw.quatern_c,'float');
                fwrite(fid,O.head.raw.quatern_d,'float');
                fwrite(fid,O.head.raw.qoffset_x,'float');
                fwrite(fid,O.head.raw.qoffset_y,'float');
                fwrite(fid,O.head.raw.qoffset_z,'float');
                fclose(fid);
            case 'Qform'
                % Add the sform
                fid = fopen(mname,'r+');
                fseek(fid,76,'bof');
                fwrite(fid,O.head.raw.pixdim(1),'float');
                fseek(fid,254,'bof');
                fwrite(fid,O.head.raw.sform_code,'short');
                fseek(fid,280,'bof');
                fwrite(fid,O.head.raw.srow_x,'float');
                fwrite(fid,O.head.raw.srow_y,'float');
                fwrite(fid,O.head.raw.srow_z,'float');
                fclose(fid);
        end
        gzip(mname)
        delete(mname)
        mname = [mname '.gz'];
    end


    function tform = matchedpoint_registration(MIP_Secondary,MIP_orig)
        % SIFT feature extraction
        featurepoints_I = detectSIFTFeatures(MIP_orig);
        [fI,vpI] = extractFeatures(MIP_orig,featurepoints_I);
        featurepoints_Secondary = detectSIFTFeatures(MIP_Secondary);
        [fS,vpS] = extractFeatures(MIP_Secondary,featurepoints_Secondary);
        % Match features
        ipA = matchFeatures(fS,fI);
        matchedpoints_Secondary = vpS(ipA(:,1));
        matchedpoints_I = vpI(ipA(:,2));
        % compute similarity transform
        try
            tform = estgeotform2d(matchedpoints_Secondary,matchedpoints_I,'similarity');
        catch
            tform = [];
            return
        end
    end


    function [cI,idx,corresponds] = find_substack(f2,f1,I)
        swmp = zeros(length(f1),length(f2));
        for s=1:length(f1)
            spoints = f1{s};
            parfor i=1:length(f2)
                ipoints = f2{i};
                ipA = matchFeatures(spoints,ipoints);
                swmp(s,i) = size(ipA,1);
            end
        end
        [mxpts,idx] = max(swmp,[],2);
        if any(mxpts==0) || ~all(abs(diff(idx))==1)
            % Check for close correspondance
            il = length(idx);
            mdl = fitlm(0:il-1,idx);
            incpt = round(mdl.Coefficients.Estimate(1));
            slope = mdl.Coefficients.Estimate(2);
            if abs((1-abs(slope)))<0.05 && mdl.Coefficients.pValue(2)<0.0001
                % If estimated slope is ~1 and the fit is strong, construct
                % new idx
                slope = round(slope);
                idx = incpt:slope:(il-1)*slope+incpt;
                try
                    cI = I(:,:,idx);
                catch
                    corresponds = false;
                    cI = [];
                    idx = [];
                    return
                end
                corresponds = true;
            else
                corresponds = false;
                cI = [];
                idx = [];
            end
        else
            cI = I(:,:,idx);
            corresponds = true;
        end
    end

    function max_corr = quick_match(V, U)
        max_corr = -Inf;
        for i = 1:(length(U) - length(V) + 1)
            subset = U(i:i+length(V)-1);
            corr_value = corr(V,subset);
            if corr_value > max_corr
                max_corr = corr_value;
            end
            % also try flipping
            corr_value = corr(flip(V),subset);
            if corr_value > max_corr
                max_corr = corr_value;
            end
        end

    end


end
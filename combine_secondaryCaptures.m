function [Secondary_I,Secondary_Mask,maskname,Masks] = combine_secondaryCaptures(maskname)


    M = cell(length(maskname),1);
    for sc=1:length(maskname)
        % Get the dicom file metadata
        mask_dcms = struct2table(fast_dcm_head(maskname{sc},{'InstanceNumber','Rows','Columns'}));
        mask_dcms = sortrows(mask_dcms,'InstanceNumber');
        % Some directories are "dirty" with combined secondary and original
        % images. Toss out any with lower pixel count
        mask_dcms(mask_dcms.Rows<max(mask_dcms.Rows),:) = [];
        % Keep only unique instance numbers (for rare cases when secondary
        % captures had been re-saved in the same directory).
        [~,ia] = unique(mask_dcms.InstanceNumber);
        mask_dcms = mask_dcms(ia,:);
        % Instance numbers should run 0-#slices
        if ~all(diff(mask_dcms.InstanceNumber)==1)
            continue
        end
        
        
        if mask_dcms.Rows(1)~=mask_dcms.Columns(1)
            continue
        end
        % Instantiate arrary
        M{sc} = uint8(zeros([mask_dcms.Rows(1) mask_dcms.Columns(1) 3 height(mask_dcms)]));
        for s=1:height(mask_dcms)
            % Get the pixel data
            M{sc}(:,:,:,s) = dicomread(mask_dcms.Filename{s});
        end
        % swap 3rd and 4th dimensions to put color channels in dim4
        M{sc} = permute(M{sc},[1 2 4 3]);
    
        % compute correlation between channels
        chancorr = corr(double(reshape(M{sc},[],3)));
        % get the index of the channel unlike the other two
        i = find(sum(chancorr>0.99999)==1);
        % ensure mask channel is 1 (red)
        if i==2
            M{sc} = M{sc}(:,:,:,[2 1 3]);
        elseif i==3
            M{sc} = M{sc}(:,:,:,[3 2 1]);
        end
    end

    nogood = cellfun(@isempty,M);
    maskname(nogood) = [];
    M(nogood) = [];

    % take only the masks with the most common dimensionality
    [r,~,~,~] = cellfun(@size,M);
    maskname(r~=mode(r)) = [];
    M(r~=mode(r)) = [];
    
    Secondary_I = mat2gray(M{1}(:,:,:,2));
    Secondary_I(M{1}(:,:,:,1)~=M{1}(:,:,:,2)) = NaN;
    channel_ratio = double(M{1}(:,:,:,1))./double(M{1}(:,:,:,2));
    Secondary_Mask = channel_ratio>=1.5 & channel_ratio<=2.5;
    % Eliminate very small holes in the mask (use this rather than imfill
    % because there might be real holes in the segmentation)
    Secondary_Mask = ~bwareaopen(~Secondary_Mask,5,6);
    % Remove any cluster that is <1% of the mask
    cc = bwconncomp(Secondary_Mask,6);
    nvox = cellfun(@length,cc.PixelIdxList);
    small_vols = (nvox./sum(nvox))<0.01;
    Secondary_Mask(vertcat(cc.PixelIdxList{small_vols})) = 0;

    
    for sc=2:length(maskname)
        im2 = mat2gray(M{sc}(:,:,:,2));
        im2(M{sc}(:,:,:,1)~=M{sc}(:,:,:,2)) = NaN;
        channel_ratio = double(M{sc}(:,:,:,1))./double(M{sc}(:,:,:,2));
        m2 = channel_ratio>=1.5 & channel_ratio<=2.5;
        m2 = ~bwareaopen(~m2,5,6);
        cc = bwconncomp(m2,6);
        nvox = cellfun(@length,cc.PixelIdxList);
        small_vols = (nvox./sum(nvox))<0.01;
        m2(vertcat(cc.PixelIdxList{small_vols})) = 0;

        sA = Secondary_I;
        sA(round(0.9*size(sA,1)):end,round(0.75*size(sA,2)):end,:) = NaN;
        mat = corr(...
            reshape(sA,[],size(sA,3)),...
            reshape(im2,[],size(im2,3)),...
            'Rows','complete');
        mat = mat>0.99;

        if ~any(mat(:))
            % no overlap, flip l/r
            im2 = flip(im2,2);
            m2 = flip(m2,2);
            mat = corr(...
                reshape(Secondary_I,[],size(Secondary_I,3)),...
                reshape(im2,[],size(im2,3)),...
                'Rows','pairwise');
            mat = mat>0.99;
            if ~any(mat(:))
                Secondary_Mask = cat(4,Secondary_Mask,false(size(Secondary_I)));
                continue
            end
        end

        % For matched slices, average together (to recover intensity in
        % masked regions
        im2inSecondary = find(any(mat));
        Secondary_slices = find(any(mat,2));
        Secondary_I(:,:,Secondary_slices) = mean(cat(4,...
            Secondary_I(:,:,Secondary_slices),...
            im2(:,:,im2inSecondary)),4,'omitmissing');
        % Add a channel to the mask
        Secondary_Mask = cat(4,Secondary_Mask,false(size(Secondary_I)));
        Secondary_Mask(:,:,Secondary_slices,sc) = m2(:,:,im2inSecondary);

        if all(~any(mat))
            % no correspondance, do not add because you don't know to
            % prepend or append or if the two stacks are contiguous
            continue
        end

        % slices in M{sc} that are not in M{1}
        toadd = find(~any(mat));
        if ~isempty(toadd)
            if isscalar(toadd) || all(diff(toadd)==1)
                % 1 slice or contiguous slices, add before or after
                if toadd(1)==1
                    % Add slices of im2 before Secondary_I
                    Secondary_I = cat(3,im2(:,:,toadd),Secondary_I);
                    % Add N slices of false before Secondary_Mask
                    Secondary_Mask = cat(3,false([size(m2,1) size(m2,2) length(toadd) sc]),Secondary_Mask);
                    % Fill those slices with the slices of m2
                    Secondary_Mask(:,:,toadd,sc) = m2(:,:,toadd);
                else
                    % Add slices of im2 after Secondary_I
                    Secondary_I = cat(3,Secondary_I,im2(:,:,toadd));
                    % Add N slices of false after Secondary_Mask
                    ss = size(Secondary_Mask,3)+1;
                    es = ss+length(toadd)-1;
                    Secondary_Mask = cat(3,Secondary_Mask,false([size(m2,1) size(m2,2) length(toadd) sc]));
                    % Fill those slices with the slices of m2
                    Secondary_Mask(:,:,ss:es,sc) = m2(:,:,toadd);
                end
            else
                % discontiguous slices, add before and after
                splitat = find(diff(toadd)~=1);
                if length(splitat)>1
                    % This should not occur, skip
                    continue
                end
                prepend = toadd(1:splitat);
                append = toadd(splitat+1:end);
                Secondary_I = cat(3,im2(:,:,prepend),Secondary_I,im2(:,:,append));
                Secondary_Mask = cat(3,false([size(m2,1) size(m2,2) length(prepend) sc]),...
                    Secondary_Mask,false([size(m2,1) size(m2,2) length(append) sc]));
                Secondary_Mask(:,:,prepend,sc) = m2(:,:,prepend);
                Secondary_Mask(:,:,append,sc) = m2(:,:,append);
            end
        end

    end

    Masks = struct();
    for n=1:length(maskname)
        Masks(n).list = maskname{n};
        Masks(n).BW = Secondary_Mask(:,:,:,n);
    end

end
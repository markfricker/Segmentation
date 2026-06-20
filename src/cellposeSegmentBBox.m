function [L, BW] = cellposeSegmentBBox(I, params)
%CELLPOSESEGMENTBBOX  Accelerated Cellpose via coarse bounding-box pre-segmentation.
%
%   [L, BW] = cellposeSegmentBBox(I, params)
%
%   Runs a fast Otsu threshold on the full frame to find object regions,
%   dilates the coarse mask to merge nearby objects, then crops each
%   resulting bounding box and passes it to cellposeSegment individually.
%   Labels are offset and stitched back into a full-frame output.
%
%   The speedup arises because Cellpose inference time scales with image
%   area, and at typical mitochondrial densities (5-15% coverage) the
%   total bbox area is much smaller than the full frame.  An empty coarse
%   mask (no objects above threshold) short-circuits Cellpose entirely.
%
%   INPUTS
%     I      : 2-D grayscale image (any numeric class).
%     params : parameter struct — same fields as cellposeSegment, plus:
%
%       .bboxMargin      (default 20 px)
%           Dilation radius applied to the coarse binary mask before
%           computing bounding boxes.  Objects within 2*bboxMargin px of
%           each other are grouped into a single crop.  Larger values
%           give more context to Cellpose but reduce the area saving.
%
%       .bboxCoarseScale (default 1.0)
%           Multiplier on the Otsu threshold for the coarse mask.
%           Values < 1 lower the threshold (more permissive, fewer missed
%           objects); values > 1 raise it (tighter, faster).  Default 1.0
%           (pure Otsu) is generally safe for mitochondria fluorescence.
%
%       .bboxMinArea     (default 5 px^2)
%           Minimum area of a coarse blob to be considered.  Smaller
%           blobs are noise and are discarded before dilation.
%
%   OUTPUTS
%     L  : uint16 label image, full frame size, labels unique across crops.
%     BW : binary mask, single precision (L > 0).
%
%   If the coarse mask is empty, both L and BW are zero-filled and
%   cellposeSegment is never called.
%
%   NOTES
%   - Requires Image Processing Toolbox (imbinarize, regionprops, imdilate).
%   - All other params fields are forwarded unchanged to cellposeSegment.
%   - Labels in L are globally unique: each crop's labels are offset by
%     the maximum label assigned so far before stitching.
%
%   See also: cellposeSegment

    % --- defaults for bbox parameters ----------------------------------------
    if nargin < 2, params = struct(); end
    if ~isfield(params, 'bboxMargin'),      params.bboxMargin      = 20;  end
    if ~isfield(params, 'bboxCoarseScale'), params.bboxCoarseScale = 1.0; end
    if ~isfield(params, 'bboxMinArea'),     params.bboxMinArea     = 5;   end

    [nY, nX] = size(I);
    L  = zeros(nY, nX, 'uint16');
    BW = zeros(nY, nX, 'single');

    % --- coarse binary mask --------------------------------------------------
    Inorm = im2single(I);
    thresh = graythresh(Inorm) * params.bboxCoarseScale;
    BWcoarse = Inorm > thresh;

    if ~any(BWcoarse(:))
        return  % empty frame — skip Cellpose entirely
    end

    % Remove tiny noise blobs before dilation.
    if params.bboxMinArea > 1
        BWcoarse = bwareaopen(BWcoarse, params.bboxMinArea);
    end

    if ~any(BWcoarse(:))
        return
    end

    % --- dilate to merge nearby objects and compute bounding boxes ----------
    se      = strel('disk', params.bboxMargin);
    BWdil   = imdilate(BWcoarse, se);
    props   = regionprops(BWdil, 'BoundingBox');

    if isempty(props)
        return
    end

    % --- run Cellpose per crop and stitch -----------------------------------
    nextLabel = uint16(1);

    for k = 1:numel(props)
        bb = props(k).BoundingBox;   % [x y w h], 1-based pixel-edge coords

        % Convert to integer pixel indices, clipped to image bounds.
        c1 = max(1,    floor(bb(1)));
        r1 = max(1,    floor(bb(2)));
        c2 = min(nX,   floor(bb(1) + bb(3) - 1) + 1);
        r2 = min(nY,   floor(bb(2) + bb(4) - 1) + 1);

        % Skip degenerate crops (can happen at image edges).
        if c2 <= c1 || r2 <= r1
            continue
        end

        Icrop = I(r1:r2, c1:c2);

        % Run Cellpose on the crop.
        try
            Lcrop = cellposeSegment(Icrop, params);
        catch ME
            warning('cellposeSegmentBBox:cropFailed', ...
                'Crop %d/%d failed: %s', k, numel(props), ME.message);
            continue
        end

        % Skip empty crops.
        if ~any(Lcrop(:))
            continue
        end

        % Offset labels so they are globally unique, then stitch.
        Lcrop(Lcrop > 0) = Lcrop(Lcrop > 0) + (nextLabel - uint16(1));
        L(r1:r2, c1:c2) = max(L(r1:r2, c1:c2), Lcrop);
        nextLabel = max(L(:)) + uint16(1);
    end

    BW = single(L > 0);
end

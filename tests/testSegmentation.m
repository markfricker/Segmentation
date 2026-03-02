classdef testSegmentation < matlab.unittest.TestCase
% testSegmentation  Unit tests for the Segmentation toolbox
%
% USAGE
%   Run all tests from the MATLAB command window:
%       results = runtests('tests/testSegmentation');
%       table(results)
%
%   Run a single test method:
%       runtests('tests/testSegmentation/testLTF_blobDetected_sauvola')
%
% COVERAGE
%   localThresholdFast  — 10 tests
%   watershedSegment    — 13 tests
%   refineSegment       — 13 tests
%
% REQUIREMENTS
%   MATLAB R2014b+ (matlab.unittest framework; activecontour for refineSegment)
%   Image Processing Toolbox (bwdist, watershed, imgaussfilt, activecontour)
%   All src/ functions on the path (added automatically by TestClassSetup).

    properties (Constant)
        Tol = single(1e-4);   % absolute tolerance for floating-point checks
    end

    % =====================================================================
    % Path setup — runs once before any test in the class
    % =====================================================================
    methods (TestClassSetup)
        function addSrcPath(tc)
            rootDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(rootDir, 'src'));
            addpath(fullfile(rootDir, 'demos'));
        end
    end

    % =====================================================================
    % Shared synthetic test images
    %
    %   uniformImage — constant mid-grey (no structure)
    %   blobImage    — centred isotropic Gaussian blob
    %   rodImage     — horizontal elongated Gaussian rod
    %   twoBlobImage — two spatially separated compact blobs (for counting)
    % =====================================================================
    methods (Static, Access = private)

        function I = uniformImage()
            I = 0.5 * ones(64, 64, 'single');
        end

        function I = blobImage()
            % Single bright Gaussian blob (sigma=5 px) at centre of 64x64.
            [xx, yy] = meshgrid(1:64, 1:64);
            I = single(exp(-((xx-32).^2 + (yy-32).^2) / (2*5^2)));
        end

        function I = rodImage()
            % Horizontal rod: Gaussian cross-section (sigma_y=3 px), 16 px long.
            [xx, yy] = meshgrid(1:64, 1:64);
            I = single(exp(-(yy-32).^2 / (2*3^2)) .* (xx >= 24 & xx <= 39));
            I = I / max(I(:));
        end

        function I = twoBlobImage()
            % Two compact blobs separated by a clear gap; should yield
            % exactly two objects after watershed segmentation.
            [xx, yy] = meshgrid(1:64, 1:64);
            b1 = exp(-((xx-20).^2 + (yy-32).^2) / (2*4^2));
            b2 = exp(-((xx-44).^2 + (yy-32).^2) / (2*4^2));
            I  = single(max(b1, b2));
        end

        function [I, L] = blobWithLabel()
            % Raw blob image + a tight initial label for refineSegment tests.
            % The label is thresholded at 0.7 (well inside the blob peak),
            % so the labelled region is smaller than the true fluorescent
            % extent.  After refinement the mask should expand toward the
            % true blob boundary.
            I  = testSegmentation.blobImage();
            BW = I > 0.7;
            L  = uint16(bwlabel(BW));
        end

    end

    % =====================================================================
    %  localThresholdFast
    % =====================================================================
    methods (Test)

        function testLTF_smoke(tc)
            % Default call must complete without error.
            BW = localThresholdFast(testSegmentation.blobImage());
            tc.verifyNotEmpty(BW);
        end

        function testLTF_outputSize(tc)
            I  = testSegmentation.blobImage();
            BW = localThresholdFast(I);
            tc.verifySize(BW, size(I));
        end

        function testLTF_outputClass(tc)
            BW = localThresholdFast(testSegmentation.blobImage());
            tc.verifyClass(BW, 'logical');
        end

        function testLTF_binaryOutput(tc)
            % Every element must be exactly 0 or 1 (logical).
            BW = localThresholdFast(testSegmentation.blobImage());
            tc.verifyTrue(all(BW(:) == 0 | BW(:) == 1));
        end

        function testLTF_thresholdSurface_size(tc)
            % T must match image size when second output is requested.
            I = testSegmentation.blobImage();
            [~, T] = localThresholdFast(I);
            tc.verifySize(T, size(I));
        end

        function testLTF_invert_flips(tc)
            % 'invert' must produce the exact complement.
            I       = testSegmentation.blobImage();
            BW_fwd  = localThresholdFast(I);
            BW_inv  = localThresholdFast(I, 'invert', true);
            tc.verifyTrue(isequal(BW_inv, ~BW_fwd));
        end

        function testLTF_blobDetected_sauvola(tc)
            % Sauvola must classify the bright blob centre as foreground.
            I  = testSegmentation.blobImage();
            BW = localThresholdFast(I, 'method', 'sauvola', ...
                                       'windowsize', 21, 'k', 0.2, 'r', 0.5);
            tc.verifyTrue(BW(32, 32));
        end

        function testLTF_allMethodsRun(tc)
            % Every named method must execute without error on a blob image.
            I       = testSegmentation.blobImage();
            methods = {'sauvola','niblack','wolf','nick','phansalkar', ...
                       'bradley','bernsen','midgrey'};
            for m = methods
                tc.verifyNotEmpty( ...
                    localThresholdFast(I, 'method', m{1}), ...
                    sprintf('Method ''%s'' returned empty output', m{1}));
            end
        end

        function testLTF_precisionDouble(tc)
            % 'double' precision option must not error and return logical BW.
            I  = testSegmentation.blobImage();
            BW = localThresholdFast(I, 'precision', 'double');
            tc.verifyClass(BW, 'logical');
        end

        function testLTF_unknownMethod_errors(tc)
            tc.verifyError( ...
                @() localThresholdFast(testSegmentation.blobImage(), ...
                                       'method', 'frobnicator'), ...
                'localThresholdFast:unknownMethod');
        end

    end

    % =====================================================================
    %  watershedSegment
    % =====================================================================
    methods (Test)

        function testWS_smoke(tc)
            % Default (marker) call must complete without error.
            [BW, L] = watershedSegment(testSegmentation.blobImage());
            tc.verifyNotEmpty(BW);
            tc.verifyNotEmpty(L);
        end

        function testWS_outputSize_BW(tc)
            I  = testSegmentation.blobImage();
            BW = watershedSegment(I);
            tc.verifySize(BW, size(I));
        end

        function testWS_outputSize_L(tc)
            I       = testSegmentation.blobImage();
            [~, L]  = watershedSegment(I);
            tc.verifySize(L, size(I));
        end

        function testWS_outputClass_BW(tc)
            BW = watershedSegment(testSegmentation.blobImage());
            tc.verifyClass(BW, 'single');
        end

        function testWS_outputClass_L(tc)
            [~, L] = watershedSegment(testSegmentation.blobImage());
            tc.verifyClass(L, 'uint16');
        end

        function testWS_BW_isBinary(tc)
            % BW must contain only 0 and 1 (no intermediate values).
            BW = watershedSegment(testSegmentation.blobImage());
            uniqueVals = unique(BW(:));
            tc.verifyTrue(all(uniqueVals == 0 | uniqueVals == 1));
        end

        function testWS_labels_nonnegative(tc)
            [~, L] = watershedSegment(testSegmentation.blobImage());
            tc.verifyGreaterThanOrEqual(min(L(:)), uint16(0));
        end

        function testWS_blobDetected_marker(tc)
            % Marker method must detect the bright blob as at least one object.
            [~, L] = watershedSegment(testSegmentation.blobImage(), ...
                                      'method', 'marker', 'threshold', 0.3);
            tc.verifyGreaterThanOrEqual(max(L(:)), uint16(1));
        end

        function testWS_twoBlobsSeparated(tc)
            % Two well-separated blobs must yield exactly two labelled objects.
            [~, L] = watershedSegment(testSegmentation.twoBlobImage(), ...
                                      'method', 'marker', ...
                                      'threshold', 0.3, ...
                                      'hMinima',   0.05, ...
                                      'minArea',   5);
            tc.verifyEqual(max(L(:)), uint16(2));
        end

        function testWS_method_distance(tc)
            % Distance-transform method must run and return binary + labels.
            [BW, L] = watershedSegment(testSegmentation.blobImage(), ...
                                       'method', 'distance', 'threshold', 0.3);
            tc.verifyClass(BW, 'single');
            tc.verifyClass(L,  'uint16');
        end

        function testWS_method_gradient(tc)
            % Gradient method must run and return binary + labels.
            [BW, L] = watershedSegment(testSegmentation.blobImage(), ...
                                       'method', 'gradient', 'threshold', 0.3);
            tc.verifyClass(BW, 'single');
            tc.verifyClass(L,  'uint16');
        end

        function testWS_minArea_removesAllObjects(tc)
            % Setting minArea larger than any object must yield empty output.
            [BW, L] = watershedSegment(testSegmentation.blobImage(), ...
                                       'threshold', 0.3, 'minArea', 10000);
            tc.verifyEqual(max(BW(:)), single(0));
            tc.verifyEqual(max(L(:)),  uint16(0));
        end

        function testWS_unknownMethod_errors(tc)
            tc.verifyError( ...
                @() watershedSegment(testSegmentation.blobImage(), ...
                                     'method', 'frobnicator'), ...
                'watershedSegment:unknownMethod');
        end

    end

    % =====================================================================
    %  refineSegment
    % =====================================================================
    methods (Test)

        function testRS_smoke_chanvese(tc)
            % Chan-Vese call must complete without error.
            [I, L] = testSegmentation.blobWithLabel();
            [BW_r, L_r] = refineSegment(I, L, 'method', 'chanvese', ...
                                               'maxExpand', 4, 'nIter', 30);
            tc.verifyNotEmpty(BW_r);
            tc.verifyNotEmpty(L_r);
        end

        function testRS_smoke_dilate(tc)
            % Dilate call must complete without error.
            [I, L] = testSegmentation.blobWithLabel();
            [BW_r, L_r] = refineSegment(I, L, 'method', 'dilate', ...
                                               'maxExpand', 4);
            tc.verifyNotEmpty(BW_r);
            tc.verifyNotEmpty(L_r);
        end

        function testRS_outputSize_BW(tc)
            % BW_r must be the same size as the input image.
            [I, L] = testSegmentation.blobWithLabel();
            BW_r = refineSegment(I, L, 'method', 'dilate');
            tc.verifySize(BW_r, size(I));
        end

        function testRS_outputSize_L(tc)
            % L_r must be the same size as the input image.
            [I, L] = testSegmentation.blobWithLabel();
            [~, L_r] = refineSegment(I, L, 'method', 'dilate');
            tc.verifySize(L_r, size(I));
        end

        function testRS_outputClass_BW(tc)
            % BW_r must be single precision.
            [I, L] = testSegmentation.blobWithLabel();
            BW_r = refineSegment(I, L, 'method', 'dilate');
            tc.verifyClass(BW_r, 'single');
        end

        function testRS_outputClass_L(tc)
            % L_r must be uint16.
            [I, L] = testSegmentation.blobWithLabel();
            [~, L_r] = refineSegment(I, L, 'method', 'dilate');
            tc.verifyClass(L_r, 'uint16');
        end

        function testRS_BW_isBinary(tc)
            % BW_r must contain only 0 and 1 (no intermediate values).
            [I, L] = testSegmentation.blobWithLabel();
            BW_r = refineSegment(I, L, 'method', 'dilate');
            uniqueVals = unique(BW_r(:));
            tc.verifyTrue(all(uniqueVals == 0 | uniqueVals == 1));
        end

        function testRS_labels_nonnegative(tc)
            % No label in L_r may be negative.
            [I, L] = testSegmentation.blobWithLabel();
            [~, L_r] = refineSegment(I, L, 'method', 'dilate');
            tc.verifyGreaterThanOrEqual(min(L_r(:)), uint16(0));
        end

        function testRS_dilate_expands(tc)
            % Dilate refinement must increase (or equal) the foreground area.
            [I, L] = testSegmentation.blobWithLabel();
            areaIn  = nnz(L > 0);
            BW_r    = refineSegment(I, L, 'method', 'dilate', ...
                                    'maxExpand', 6, 'fgThreshold', 0.05);
            areaOut = nnz(BW_r);
            tc.verifyGreaterThanOrEqual(areaOut, areaIn);
        end

        function testRS_chanvese_expands(tc)
            % Chan-Vese refinement must increase (or equal) the foreground area.
            [I, L] = testSegmentation.blobWithLabel();
            areaIn  = nnz(L > 0);
            BW_r    = refineSegment(I, L, 'method', 'chanvese', ...
                                    'maxExpand', 6, 'nIter', 50);
            areaOut = nnz(BW_r);
            tc.verifyGreaterThanOrEqual(areaOut, areaIn);
        end

        function testRS_dilate_maxExpand0_noGrowth(tc)
            % dilate with maxExpand=0 must not add any new foreground pixels
            % (only original seeds are valid; no expansion zone exists).
            [I, L] = testSegmentation.blobWithLabel();
            areaIn  = nnz(L > 0);
            BW_r    = refineSegment(I, L, 'method', 'dilate', 'maxExpand', 0);
            areaOut = nnz(BW_r);
            tc.verifyLessThanOrEqual(areaOut, areaIn);
        end

        function testRS_emptyLabel_returnsZeros(tc)
            % An all-zero label image must return all-zero BW_r and L_r.
            I   = testSegmentation.blobImage();
            L   = zeros(size(I), 'uint16');
            [BW_r, L_r] = refineSegment(I, L, 'method', 'dilate');
            tc.verifyEqual(max(BW_r(:)), single(0));
            tc.verifyEqual(max(L_r(:)),  uint16(0));
        end

        function testRS_preservesObjectCount_dilate(tc)
            % Dilate refinement must not increase the number of labelled objects
            % (boundaries expand; no new objects are created or merged).
            [I, L] = testSegmentation.blobWithLabel();
            nIn     = double(max(L(:)));
            [~, L_r] = refineSegment(I, L, 'method', 'dilate', 'maxExpand', 4);
            nOut    = double(max(L_r(:)));
            tc.verifyEqual(nOut, nIn);
        end

        function testRS_unknownMethod_errors(tc)
            [I, L] = testSegmentation.blobWithLabel();
            tc.verifyError( ...
                @() refineSegment(I, L, 'method', 'frobnicator'), ...
                'refineSegment:unknownMethod');
        end

    end

end

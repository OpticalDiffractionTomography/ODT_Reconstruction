classdef tHelpers < matlab.unittest.TestCase
    % Unit tests for the numerical and image helpers in src/:
    % mk_ellipse, PhiShift, phaseCompensation, saveTomogramTIFF,
    % saveTomogramPNG, saveInspectionPNG.
    %   results = runtests('tests')   (run from the repo root)
    % Needs the Image Processing Toolbox (used by the PNG/TIFF helpers).

    properties
        Dir
    end

    methods (TestClassSetup)
        function addSrcToPath(testCase)
            srcDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(srcDir));
        end
    end

    methods (TestMethodSetup)
        function makeTempDir(testCase)
            testCase.Dir = tempname;
            mkdir(testCase.Dir);
            testCase.addTeardown(@() rmdir(testCase.Dir, 's'));
        end
    end

    methods (Test)
        % ── mk_ellipse ────────────────────────────────────────────────────
        function ellipseSizeIsYByX(testCase)
            H = mk_ellipse(10, 10, 80, 60);
            testCase.verifySize(H, [60 80]);
            testCase.verifyClass(H, 'logical');
        end

        function ellipseIsFalseInsideAndTrueOutside(testCase)
            % true marks pixels OUTSIDE the ellipse; the centre is (X/2, Y/2)
            H = mk_ellipse(10, 10, 64, 64);
            testCase.verifyFalse(H(32, 32));      % centre
            testCase.verifyFalse(H(32, 42));      % exactly on the radius: inside
            testCase.verifyTrue(H(32, 43));       % one pixel beyond
            testCase.verifyTrue(H(1, 1));         % corner
        end

        function ellipseRespectsDifferentRadii(testCase)
            H = mk_ellipse(20, 5, 100, 100);
            testCase.verifyFalse(H(50, 50 + 20));   % along x, within XR
            testCase.verifyTrue(H(50 + 6, 50));     % along y, beyond YR
        end

        % ── PhiShift ──────────────────────────────────────────────────────
        function phiShiftLeavesFlatImageUnchanged(testCase)
            U = 5 * ones(40, 50);
            testCase.verifyEqual(PhiShift(U), U, 'AbsTol', 1e-12);
        end

        function phiShiftRemovesMostOfALinearTilt(testCase)
            [X, Y] = meshgrid(1:200, 1:200);
            U = 0.05 * X - 0.03 * Y + 2;
            Out = PhiShift(U);
            % the border strips are 4 px wide, so a small residual remains
            testCase.verifyLessThan(range(Out(:)), 0.1 * range(U(:)));
        end

        function phiShiftKeepsSize(testCase)
            testCase.verifySize(PhiShift(rand(30, 45)), [30 45]);
        end

        % ── phaseCompensation ────────────────────────────────────────────
        function degreeOneRemovesAPlane(testCase)
            [X, Y] = meshgrid(1:20, 1:20);
            p = 50 + 0.1 * X - 0.2 * Y;
            [good, coef, comp] = phaseCompensation(p, 1);
            testCase.verifyEqual(good, zeros(20), 'AbsTol', 1e-6);
            testCase.verifyNumElements(coef, 3);        % 2n + 1
            testCase.verifySize(comp, [20 20]);
        end

        function degreeTwoRemovesAQuadraticSurface(testCase)
            [X, Y] = meshgrid(1:20, 1:20);
            p = 50 + 0.01 * X.^2 + 0.02 * Y.^2 + 0.003 * X - 0.004 * Y;
            [good, coef] = phaseCompensation(p, 2);
            testCase.verifyEqual(good, zeros(20), 'AbsTol', 1e-4);
            testCase.verifyNumElements(coef, 5);        % 2n + 1
        end

        function goodPlusCompensationRecoversInput(testCase)
            p = 10 + rand(20);
            [good, ~, comp] = phaseCompensation(p, 2);
            testCase.verifyEqual(good + comp, p, 'AbsTol', 1e-10);
        end

        function maskExcludesCellsFromTheFit(testCase)
            [X, Y] = meshgrid(1:20, 1:20);
            bg = 50 + 0.1 * X - 0.2 * Y;
            p = bg;  p(8:12, 8:12) = p(8:12, 8:12) + 10;       % a "cell"
            mask = true(20);  mask(8:12, 8:12) = false;
            outside = ~ismember(1:20, 8:12);

            masked = phaseCompensation(p, 1, mask);
            testCase.verifyEqual(masked(outside, outside), zeros(15), 'AbsTol', 1e-4);
            testCase.verifyEqual(masked(8:12, 8:12), 10 * ones(5), 'AbsTol', 1e-4);

            unmasked = phaseCompensation(p, 1);     % the cell pulls the fit off
            testCase.verifyGreaterThan(max(abs(unmasked(outside, outside)), [], 'all'), 0.05);
        end

        % ── saveTomogramTIFF ─────────────────────────────────────────────
        function tiffHasOnePagePerSliceAndFlipsLeftRight(testCase)
            R = 1.33 + 0.07 * rand(8, 6, 3);
            f = fullfile(testCase.Dir, 'tomo.tif');
            saveTomogramTIFF(R, f);
            testCase.verifyNumElements(imfinfo(f), 3);
            for k = 1:3
                expected = fliplr(uint16(real(R(:, :, k)) * 10000));
                testCase.verifyEqual(imread(f, k), expected);
            end
        end

        function tiffIsOverwrittenNotAppendedOnRerun(testCase)
            R = 1.33 + 0.07 * rand(8, 6, 3);
            f = fullfile(testCase.Dir, 'tomo.tif');
            saveTomogramTIFF(R, f);
            saveTomogramTIFF(R, f);
            testCase.verifyNumElements(imfinfo(f), 3);
        end

        function tiffIgnoresImaginaryPart(testCase)
            R = (1.33 + 0.07 * rand(8, 6, 2)) + 1i * rand(8, 6, 2);
            f = fullfile(testCase.Dir, 'tomo.tif');
            saveTomogramTIFF(R, f);
            testCase.verifyEqual(imread(f, 1), fliplr(uint16(real(R(:, :, 1)) * 10000)));
        end

        % ── PNG helpers ──────────────────────────────────────────────────
        function tomogramPngIsWritten(testCase)
            R = 1.337 + 0.03 * rand(32, 32, 16);
            f = fullfile(testCase.Dir, 'tomo.png');
            saveTomogramPNG(R, 1.377, f);
            testCase.verifyTrue(isfile(f));
            testCase.verifyGreaterThan(numel(imread(f)), 0);
        end

        function inspectionPngIsWrittenForFieldFiles(testCase)
            retPhase = 0.1 * rand(10, 10, 5); %#ok<NASGU>
            save(fullfile(testCase.Dir, 'Field_a.mat'), 'retPhase');
            save(fullfile(testCase.Dir, 'Field_b.mat'), 'retPhase');
            f = fullfile(testCase.Dir, 'Field_inspection.png');
            saveInspectionPNG(testCase.Dir, dir(fullfile(testCase.Dir, 'Field*.mat')), f);
            testCase.verifyTrue(isfile(f));
            testCase.verifySize(imread(f), [400 800 3]);
        end

        function inspectionPngToleratesAnEmptyFieldList(testCase)
            f = fullfile(testCase.Dir, 'Field_inspection.png');
            saveInspectionPNG(testCase.Dir, dir(fullfile(testCase.Dir, 'Field*.mat')), f);
            testCase.verifyTrue(isfile(f));
        end
    end
end

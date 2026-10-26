classdef tLoadMatSafe < matlab.unittest.TestCase
    % Unit tests for loadMatSafe and recordSkipped.
    %   results = runtests('tests')   (run from the repo root)

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
        function loadsValidFile(testCase)
            f = fullfile(testCase.Dir, 'ok.mat');
            tomogMap = uint8(magic(4)); NA = 1; %#ok<NASGU>
            save(f, 'tomogMap', 'NA');
            [S, err] = loadMatSafe(f, {'tomogMap', 'NA'});
            testCase.verifyEmpty(err);
            testCase.verifyEqual(S.tomogMap, uint8(magic(4)));
        end

        function reportsMissingFile(testCase)
            [S, err] = loadMatSafe(fullfile(testCase.Dir, 'nope.mat'), {'tomogMap'});
            testCase.verifyEmpty(S);
            testCase.verifyNotEmpty(err);
        end

        function reportsGarbageFile(testCase)
            f = fullfile(testCase.Dir, 'garbage.mat');
            fid = fopen(f, 'w'); fwrite(fid, uint8(randi(255, 1, 500))); fclose(fid);
            [S, err] = loadMatSafe(f, {'tomogMap'});
            testCase.verifyEmpty(S);
            testCase.verifyNotEmpty(err);
        end

        function reportsTruncatedFile(testCase)
            % The real-world case: a copy that stopped part-way through.
            f = fullfile(testCase.Dir, 'trunc.mat');
            tomogMap = rand(200, 200, 10); %#ok<NASGU>
            save(f, 'tomogMap');
            bytes = fileread_bytes(f);
            fid = fopen(f, 'w'); fwrite(fid, bytes(1:floor(end/2))); fclose(fid);
            [S, err] = loadMatSafe(f, {'tomogMap'});
            testCase.verifyEmpty(S);
            testCase.verifyNotEmpty(err);
        end

        function reportsMissingVariables(testCase)
            f = fullfile(testCase.Dir, 'partial.mat');
            tomogMap = 1; %#ok<NASGU>
            save(f, 'tomogMap');
            [S, err] = loadMatSafe(f, {'tomogMap', 'NA', 'lambda'});
            testCase.verifyEmpty(S);
            testCase.verifySubstring(err, 'NA');
            testCase.verifySubstring(err, 'lambda');
            testCase.verifyEmpty(strfind(err, 'tomogMap'));  % present, so not listed
        end

        function neverThrows(testCase)
            testCase.verifyWarningFree(@() loadMatSafe(testCase.Dir, {'x'}));  % a folder, not a file
        end

        function recordSkippedAppendsTabSeparatedLines(testCase)
            recordSkipped(testCase.Dir, 'a.mat', 'field_Retrieval', 'cannot read');
            recordSkipped(testCase.Dir, 'b.mat', 'tomogram_Reconstruction', 'processing failed');
            lines = splitlines(strtrim(fileread(fullfile(testCase.Dir, 'skipped_files.txt'))));
            testCase.verifyNumElements(lines, 2);
            testCase.verifyEqual(strsplit(lines{1}, "\t"), {'a.mat', 'field_Retrieval', 'cannot read'});
            testCase.verifyNumElements(strsplit(lines{2}, "\t"), 3);
        end

        function recordSkippedStripsNewlinesFromReason(testCase)
            recordSkipped(testCase.Dir, 'a.mat', 'stage', sprintf('multi\nline'));
            lines = splitlines(strtrim(fileread(fullfile(testCase.Dir, 'skipped_files.txt'))));
            testCase.verifyNumElements(lines, 1);
        end

        function recordSkippedToleratesUnwritableDir(testCase)
            testCase.verifyWarningFree(@() recordSkipped(fullfile(testCase.Dir, 'missing'), 'a', 's', 'r'));
        end
    end
end

function bytes = fileread_bytes(f)
    fid = fopen(f, 'r'); bytes = fread(fid, Inf, '*uint8'); fclose(fid);
end

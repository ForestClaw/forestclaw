    filename = "fort_frame_0000.vtu";
    fid = fopen(filename, 'r');

    % Search for NumberOfPoints and NumberOfCells
    content = fscanf(fid, '%c', inf);
    numPointsStr = regexp(content, 'NumberOfPoints="(\d+)"', 'tokens');
    numCellsStr = regexp(content, 'NumberOfCells="(\d+)"', 'tokens');

    numPoints = str2double(numPointsStr{1}{1});
    numCells = str2double(numCellsStr{1}{1});

    % Search for offsets
    offsetStr = regexp(content, 'offset="(\d+)"', 'tokens');
    offsets = cellfun(@str2double, [offsetStr{:}]);

    fclose(fid);
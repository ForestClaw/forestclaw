classdef TextParser < handle
    properties
        lines
        curr_pos
    end
    
    methods
        function obj = TextParser(lines)
            obj.lines = lines;
            obj.curr_pos = 1;
        end
        
        function tokens = next(obj)
            tokens = strsplit(strtrim(obj.lines{obj.curr_pos}), ' ');
            obj.curr_pos = obj.curr_pos + 1;
        end
        
        function tf = hasNext(obj)
            tf = obj.curr_pos <= length(obj.lines);
        end
        
        function tf = currIsEmpty(obj)
            tf = isempty(strtrim(obj.lines{obj.curr_pos}));
        end
    end
end

classdef ReadClawAscii
    methods (Static)
        function patch = parsePatch(parser)
            patch = Patch();
            patch.grid_number = str2double(parser.next{1});
            patch.amr_level = str2double(parser.next{1});
            patch.block_number = str2double(parser.next{1});
            
            % Skip mpi rank
            parser.next();
            
            patch.mx = str2double(parser.next{1});
            patch.my = str2double(parser.next{1});
            patch.xlow = str2double(parser.next{1});
            patch.ylow = str2double(parser.next{1});
            patch.dx = str2double(parser.next{1});
            patch.dy = str2double(parser.next{1});
            
            % Initialize array
            patch.q = zeros(patch.mx * patch.my, 1);
            
            % Read in data
            for i = 1:patch.mx * patch.my
                while parser.currIsEmpty()
                    parser.next();
                end
                patch.q(i) = str2double(parser.next{1});
            end
        end
        
        function readFile(file, callback)
            disp(file.name);
            patches = Patch.empty();
            fileID = fopen(file.name, 'r');
            lines = textscan(fileID, '%s', 'Delimiter', '\n');
            fclose(fileID);
            lines = lines{1};
            
            parser = TextParser(lines);
            while parser.hasNext()
                if parser.currIsEmpty()
                    parser.next();
                else
                    patch = ReadClawAscii.parsePatch(parser);
                    patches(end + 1) = patch;
                end
            end
            callback(patches);
        end
    end
end

classdef Patch
    properties
        grid_number
        amr_level
        block_number
        mx
        my
        xlow
        ylow
        dx
        dy
        q
    end
end

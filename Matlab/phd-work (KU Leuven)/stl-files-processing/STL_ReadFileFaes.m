function [F, V, fpath, fpos] = STL_ReadFile(fpath, cleanup, title, read, skip, pos)
%STL_ReadFile  Read STL file.
%
%   Syntax:
%    [F, V, fpath, fpos] = ...
%        STL_ReadFile(fpath, cleanup, title, read, skip, pos)
%
%   Input:
%    fpath:   String containing a path to the file that is to be read. If
%             set to an empty string, the path will be obtained through a
%             file input dialog. Optional, defaults to an empty string.
%    cleanup: Logical indicating whether or not redundancy should be
%             removed from the mesh. Optional, defaults to true.
%    title:   Window title for the file input dialog. Only used when fpath
%             is set to an empty string. Optional, defaults to
%             'Open STL file' if omitted or empty.
%    read:    Maximum number of triangles to read from the file. Optional,
%             defaults to Inf.
%    skip:    Number of triangles or bytes to skip before reading.
%             Optional, defaults to 0.
%    pos:     Logical indicating whether skip is a number of bytes (true)
%             or a number of triangles (false). In case of ASCII STL files,
%             it is more efficient to specify a number of bytes. Optional,
%             defaults to true.
%
%   Output:
%    F:     M-by-3 array defining a triangle mesh. The rows correspond to
%           different triangles and the columns correspond to the three
%           vertices that make up each triangle. The elements are row
%           indices into V. Will be empty if the user closes the file open
%           dialog box.
%    V:     N-by-3 array defining vertices. The rows correspond to
%           different vertices and the columns correspond to X-, Y- and Z-
%           coordinates. The elements are coordinate values. Will be empty
%           if the user closes the file open dialog box.
%    fpath: String containing a path to the file that was read. Will be
%           empty if the user closes the file open dialog box.
%    fpos:  File position indicator after reading. May be used as skip
%           input argument in bytes. Will be set to 0 if the user closes
%           the file open dialog box.
%
%
%   Effect: This function will read a triangulated mesh from an STL file.
%   It automatically detects whether the specified file is an ASCII or
%   binary file. The read, skip and pos arguments can be used to read a
%   "chunk" from an STL file, which is useful for very large files that
%   would otherwise cause memory problems.
%
%   Example:
%    %Read an STL file in chunks of 1000 triangles
%    [F{1}, V{1}, fpath, fpos] = STL_ReadFile('', true, 'STL', 1000);
%    while ~isempty(F{end})
%        [F{end+1,1}, V{end+1,1}, fpath, fpos] = ...
%            STL_ReadFile(fpath, true, '', 1000, fpos, true);
%    end
%    [F, V] = TRI_Merge(F(1:end-1), V(1:end-1), true);
%
%   Dependencies: TRI_RemoveInvalidTriangles.m
%                 TRI_RemoveBadlyConnectedTriangles.m
%
%   Known parents: TRI_RegionCutter.m
%                  TRI_CutWithPlane.m
%                  TRI_CutWithMultiPlane.m
%                  Muscle_DatabaseLandmarks.m
%                  Muscle_DatabaseAllRegions.m
%                  Muscle_SelectDatabaseBone.m
%                  Muscle_ProjectDatabaseRegions.m

%Created on 10/08/2006 by Ward Bartels.
%WB, 06/02/2007: Added customizable window title.
%WB, 15/02/2007: Added badly connected triangle check and chunked reading.
%Stabile, fully functional.


%---------------%
% Process input %
%---------------%

%Set defaults for input variables
if nargin<1, fpath = ''; end
if nargin<2, cleanup = true; end
if nargin<3 || isempty(title), title = 'Open STL file'; end
if nargin<4, read = Inf; end
if nargin<5, skip = 0; end
if nargin<6, pos = true; end

%If no file path was provided, get file from file input dialog
if isempty(fpath)
    [filename, pathname] = uigetfile({'*.stl' 'STL files (*.stl)'}, title);
    if isequal(filename, 0)
        F = zeros(0, 3);
        V = zeros(0, 3);
        fpos = 0;
        return
    end
    fpath = fullfile(pathname, filename);
end


%-----------%
% Read file %
%-----------%

%Catch errors while reading, attempt to close file if one was caught
try
    
%     %Open file
%     fid = fopen(fpath, 'r');
%     
%     %Read file header and look for ASCII signature
%     header = fread(fid, 84, '*char').';
%     loc = regexpi(header, 'solid.*(facet)\s*normal', 'tokenExtents');
    
    %Check if header looks like ASCII file
    if ~isempty(loc)
        
%         %Skip specified number of bytes or triangles
%         if pos && skip~=0
%             fseek(fid, skip, 'bof');
%         else
%             fseek(fid, loc{1}(1)-1, 'bof');
%             textscan(fid, '%*s%*s%*f%*f%*f%*s%*s%*s%*f%*f%*f%*s%*f%*f%*f%*s%*f%*f%*f%*s%*s', ...
%                      skip, 'whitespace', ' \b\t\n\r');
%         end
        
        %Handle Inf value for read
        if isfinite(read)
            read = {read};
        else
            read = {};
        end
        
        %Read data from file for specified number of triangles
        [V, fpos] = textscan(fid, '%*s%*s%*f%*f%*f%*s%*s%*s%f%f%f%*s%f%f%f%*s%f%f%f%*s%*s', ...
                             read{:}, 'whitespace', ' \b\t\n\r');
        
        %Close file
        fclose(fid);
        
        %Extract vertices from data
        V = reshape([V{:}].', 3, []).';
        
    else %Read binary
        
        %Skip specified number of bytes or triangles
        if ~pos
            skip = 50*skip;
        end
        fseek(fid, skip, 'cof');
        
        %Read data from file for specified number of triangles
        V = fread(fid, [12 read], '12*float32', 2);
        
        %Save file position indicator and close file
        fpos = ftell(fid)-84;
        fclose(fid);
        
        %Extract vertices from data
        V = reshape(V(4:end,:), 3, []).';
    end
    
catch
    err = lasterror;
    try
        fclose(fid);
    end
    rethrow(err);
end


%-------------------%
% Remove redundancy %
%-------------------%

%Remove redundancy if requested, otherwise create trivial F
if cleanup
    
    %Remove redundant vertices
    [V, ignoble, jnd] = unique(V, 'rows');
    F = reshape(jnd, 3, []).';
    
    %Remove redundant triangles
    [ignoble, ind] = unique(sort(F, 2), 'rows');
    F = F(ind,:);
    
    %Remove invalid triangles <<TRI_RemoveInvalidTriangles.m>>
    F = TRI_RemoveInvalidTriangles(F);
    
    %Remove badly connected triangles and unused vertices
    %<<TRI_RemoveBadlyConnectedTriangles.m>>
    [F, V] = TRI_RemoveBadlyConnectedTriangles(F, V);
    
else
    F = reshape(1:size(V, 1), 3, []).';
end
function gitInfo = getGitInfo(gitdir)
% GITINFO = GETGITINFO(MAINDIR) - try to get information about the Git repository in MAINDIR
% (or pwd/the same directory as this script), including: 
%     - branch name of the current Git Repo 
%     - Git SHA1 HASH of the most recent commit
%
% The function first checks to see if a .git/ directory is present. If so it
% reads the .git/HEAD file to identify the branch name and then it looks up
% the corresponding commit.
%
% It then reads the .git/config file to find out the url of the
% corresponding remote repository. This is all stored in a gitInfo struct.
%
% Note this uses only file information, it makes no external program 
% calls at all. 
%
% This function must be in the base directory of the git repository
%
% Released under a BSD open source license. Based on a concept by Marc
% Gershow.
%
% Andrew Leifer
% Harvard University
% Program in Biophysics, Center for Brain Science, 
% and Department of Physics
% leifer@fas.harvard.edu
% http://www.andrewleifer.com
% 12 September 2011 
%
% Copyright 2011 Andrew Leifer. All rights reserved.
% MODIFIED: MHA 04/2018 Capdevila ITE
% 
% Redistribution and use in source and binary forms, with or without modification, are
% permitted provided that the following conditions are met:
% 
%    1. Redistributions of source code must retain the above copyright notice, this list of
%       conditions and the following disclaimer.
% 
%    2. Redistributions in binary form must reproduce the above copyright notice, this list
%       of conditions and the following disclaimer in the documentation and/or other materials
%       provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY <COPYRIGHT HOLDER> ''AS IS'' AND ANY EXPRESS OR IMPLIED
% WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
% FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
% ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% The views and conclusions contained in the software and documentation are those of the
% authors and should not be interpreted as representing official policies, either expressed
% or implied, of <copyright holder>.

headfile = '.git/HEAD';

gitInfo=struct('branch',NaN,'hash',NaN);

if nargin > 0
    switch exist(gitdir,'file')
        case 2, headfile = gitdir; % already a file, do nothing
        case 7
            if exist(fullfile(gitdir,headfile),'file'), headfile = fullfile(gitdir,headfile);
            elseif exist(fullfile(gitdir,'HEAD'),'file'), headfile = fullfile(gitdir,'HEAD');
            else, headfile = '';
            end
        otherwise, headfile = '';
    end
else
% Search up to two directories above of pwd and this script's directory
    defaults = {fileparts(mfilename('fullpath'));pwd()};
    defaults = [{defaults};cellfun(@(x) absolutepath({'..';'../..'},x),defaults,'unif',0)];
    defaults = unique(cat(2,defaults{:}),'stable');
    defaults = cellfun(@(d) fullfile(d,headfile),defaults,'unif',0);
    headfile = '';
    for j = 1:numel(defaults)
        if exist(defaults{j},'file'), headfile = defaults{j}; break; end
    end
end

if isempty(headfile)
    warning('getGitInfo:nofile','Failed to find .git/HEAD file');
    return;
end

%Read in the HEAD information, this will tell us the location of the file
%containing the SHA1
text=fileread(headfile);
parsed=textscan(text,'%s');

if ~strcmp(parsed{1}{1},'ref:') || ~length(parsed{1})>1
    warning('getGitInfo:badfile','Unexpected file format (%s)',relativepath(headfile));
    return;
end

path=parsed{1}{2};
[pathstr, name, ext]=fileparts(path);
branchName=name;

%save branchname
gitInfo.branch=branchName;

%Read in SHA1
SHA1text=fileread(fullfile(fileparts(headfile),pathstr,name,ext));
SHA1=textscan(SHA1text,'%s');
gitInfo.hash=SHA1{1}{1};
function main2

clear all; close all;


%========================================================
%Starting conformation
fname = '../common2/starting.pdb';

%Target Conformation
fnamec = '../target.pdb';

% Read those pdb data

atoms = readPDBfull(fname,[],0);   
atomsc = readPDBfull(fnamec,[],0);   

% Find the C-alphas

caidx = find(strcmp(atoms.name,'CA'));
caidxc = find(strcmp(atomsc.name,'CA'));

%========================================================
% Parameters for ANM

afftyPar.minDist = 0;  % angstroms 
%Cut-off distance
afftyPar.maxDist = 15.0; % angstroms  
afftyPar.debug   = 0;
afftyPar.dim     = 3;

%Generate the Hessian matrices
[Hc,coords] = hessCalpha_mg(atoms,afftyPar,caidx);
[Hcc,coordsc] = hessCalpha_mg(atomsc,afftyPar,caidxc);

%Coords is the coordinate matrix in the following form:
%x 1 2 3 4
%y 1 2 3 4
%z 1 2 3 4

%Perform single value decomposition to obtain the ANM modes

[Uc,Sc,Vc] = svd(full(Hc));
Sc  = diag(Sc);

% Note that the last 6 eigenvalues are zero
%Therefore when taking their inverse the last six eigenvalues are discarded

e=1./Sc(1:end-6,:).^0.5;


%Set the probabilities of each eigenvector by appending them intervals on a scale from 0 to 1

et=sum(e);
e2=e;
e2=e2/et;

N(1)=0;
for b=1:size(e2,1)
    N(b+1)=N(b)+e2(b);
end

% Initial coordinates

coordsini =coords;


%Number of CA atoms
siz1=size(coords,2);

%Evaluate  Deviation
A4=coordsc-coords;


%Select stepsize for all MC/Metropolis steps combined.
%0.5 is how much, in terms of RMSD, the end conformer will deviate from the 
%starting conformer. i.e the total deviation in the MC/Metropolis scheme

stepsiz=0.5*(siz1^0.5);

%Scaling parameter for the step in terms of square deviation
devi=0.1;

%Step size is devi scaled by the eigenvalue of the first ANM mode.
S=abs(devi*Sc(end-6))^0.5;

say=0;
say2=0;
say3=0;

%MC condition
prevper=load('oran.dat');

% a is the current parameter inside the MC/Metropolis condition
a=prevper(5);

% Here %90 is the aimed MC/Metropolis acceptance ratio
% Everything that falls outside the %85-95 range needs to be readjusted as
% follows

if prevper(2)>0.95
    a=prevper(5)*1.5;
elseif prevper(2)<0.85
    a=prevper(5)/1.5;
else
    fid2=fopen('../status.dat','w');
    fprintf(fid2,'%d',1);
    fclose(fid2);

end

%Now the initial inter-residue distances are calculated
nc=zeros(siz1,siz1);
for i1=1:siz1
    for i2=i1+1:siz1
        nc(i1,i2)=norm(coordsc(:,i1)-coordsc(:,i2) );
    end
end

%Evaluate the initial deviation energy from the target
Ep=zeros(siz1,siz1);
for i1=1:siz1
    for i2=i1+1:siz1
        Ep=Ep+(norm(coords(:,i1)-coords(:,i2)) - nc(i1,i2))^2;
    end
end

%Now we start performing the MC/Metropolis steps
%Please note the maximum step number is 100000

for b=1:100000
    b
    
    %Select a mode probabilistically
    rn=rand();
    [ind]=find(N<rn);
    ID=ind(end);
   

    %Initialize the new coordinate matrix
    coordsnew=coords;
    
    %There are 2 directions +/-
    
    %Add eigenvevtors to each row of the coordinate matrix
    %Select direction randomly
    
    %Preallocate variable for speed
    coordsnewb=zeros(3,siz1);
    
    if randi([0 1])==1
        for j=1:3
            coordsnewb(j,:)=coordsnew(j,:)+Uc(j:3:end,ID)'*e(ID)*S;
        end        
    else
        for j=1:3
            coordsnewb(j,:)=coordsnew(j,:)-Uc(j:3:end,ID)'*e(ID)*S;
        end
    end
    
    %Calculate new deviation energy from target
   
    En=0;
    for i1=1:siz1
        for i2=i1+1:siz1
            En=En+(norm(coordsnewb(:,i1)-coordsnewb(:,i2)) - nc(i1,i2))^2;
        end
    end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   %%%%%%%%%%%%MONTE CARLO CONDITION%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    if Ep>En
        say3=say3+1;
        coordsnew=coordsnewb;            
       % Accept   
       Ep=En; 
       
       
    elseif exp(-(En-Ep)*a)>rand()
        coordsnew=coordsnewb;
        say=say+1;
        say2=say2+1;
        %Accept uphill movement
        Ep=En;
    else
        %Reject step
        say=say+1;  
    end
    
    
    coords=coordsnew;
    
    
    %Check convergence conditions
    
    A3=coordsnew-coordsini;

    if sum(A3(:).^2)^0.5 > stepsiz
        break
    end
    
    X(b,:)=coordsnew(1,:);
    Y(b,:)=coordsnew(2,:);
    Z(b,:)=coordsnew(3,:);

end



writedcd('final_structure.dcd',coordsnew(1,:)' ,coordsnew(2,:)',coordsnew(3,:)');

writedcd('cg.dcd',X' ,Y',Z');


b

say2/say

say2/b


fid2=fopen('oran.dat','w');
fprintf(fid2,'%f %f %f %d %f',say2/b,say2/say,say2,b,a);
fclose(fid2);

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
function ATOM = readPDBfull(fname1,Fields,OrigPDB)
%dbstop if error
%this function d=reads PDB file into structure array
%input:
%fname1 - filename for PDB file
%Fields	- Vector indicating fields of PDB to read (default - all)
%KeepPDB - Keep full PDB content
%Output:
%ATOM - structure array with the fields corresponding to ones in ATOM record of PDB file


%TempZipFileName=unzip_AGA(fname1);  %Temporary Unzip file, if necesary

switch nargin
 case 1
  OrigPDB=0;	%Don't use information from original file
  Fields=ones(1,15);
 case 2
  OrigPDB = 0;
  if length(Fields)==0; Fields=ones(1,15); end
 case 3
  if length(Fields)==0; Fields=ones(1,15); end
end


%% the way the code was
%if nargin == 1
%	Fields=ones(1,15);
%	OrigPDB=1;
%end

%fid = fopen (fname1,'r');
%numlines = 0; 
%ENDT = ' ';
%====ATOM strings reading
%while ENDT ~= 'END'
%	numlines = numlines + 1;
%	t = fgets (fid);
%	A (numlines,1:80)=' ';
%	A (numlines, 1:length(t)) = t;
%	ENDT = A (numlines, 1:3);
%end 
%fclose (fid);

A=textread(fname1,'%s','delimiter','\n','whitespace','');
numlines=length(A);
t(1:80)=' '; t(1:length(A{1}))=A{1};
A{1}=t;
A=char(A);

%============ATOM array initialization
ATOM.type={};   % type: atom or hetatom
ATOM.line=[];	%Number of the line in PDB file
ATOM.index=[];	%        Atom serial number.
ATOM.name={};	%          Atom name.
ATOM.altloc={};	%        Alternate location indicator.
ATOM.resname={};	%       Residue name.
ATOM.chain={};	%       Chain identifier.
ATOM.resid=[];	%        Residue sequence number.
ATOM.icode={};	%         Code for insertion of residues.
ATOM.x=[];	%             Orthogonal coordinates for X in Angstroms.
ATOM.y=[];	%             Orthogonal coordinates for Y in Angstroms.
ATOM.z=[];	%             Orthogonal coordinates for Z in Angstroms.
ATOM.occupancy=[];	%     Occupancy.
ATOM.beta=[];	%    Temperature factor.
ATOM.segname={};	%         Segment identifier, left-justified.
ATOM.element={};	%       Element symbol, right-justified.
ATOM.charge={};	%        Charge on the atom.

if OrigPDB==1	%Remember the content of original PDB file
	ATOM.pdb=A;	%full content of the PDB file
else
	ATOM.pdb=[];	%no content of the PDB file
end

%  1 -  6        Record name     &quot;ATOM  &quot;
%  7 - 11        Integer         serial        Atom serial number. (1)
% 13 - 16        Atom            name          Atom name. (2)
% 17             Character       altLoc        Alternate location indicator.(3)
% 18 - 20        Residue name    resName       Residue name. (4)
% 22             Character       chainID       Chain identifier. (5)
% 23 - 26        Integer         resSeq        Residue sequence number. (6)
% 27             AChar           iCode         Code for insertion of residues.(7)
% 31 - 38        Real(8.3)       x             Orthogonal coordinates for X in
%                                              Angstroms. (8)
% 39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in
%                                              Angstroms. (9)
% 47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in
%                                              Angstroms. (10)
% 55 - 60        Real(6.2)       occupancy     Occupancy. (11)
% 61 - 66        Real(6.2)       tempFactor    Temperature factor. (12)
% 73 - 76        LString(4)      segID         Segment identifier, left-justified.
% 77 - 78        LString(2)      element       Element symbol, right-justified.
% 79 - 80        LString(2)      charge        Charge on the atom.

tic; ProgressReportTime=toc;

%=========== ATOM array filling out
for i = 1 : numlines
  %if (A(i,1:4)=='ATOM' | A(i,1:4)=='HETA')
  if (A(i,1:4)=='ATOM')    
    ATOM.line(end+1)=i;
    % Atom type: atom or heta    
    ATOM.type{end+1}=sscanf(A(i,1:4),'%s');			      
    % Atom serial number      
    if Fields(1)==1; ATOM.index(end+1)=sscanf(A(i,7:11),'%i'); end	
							
    % Atom name.
    if Fields(2)==1; ATOM.name{end+1}=sscanf(A(i,13:16),'%s'); end	
    
    % Alternate location indicator.
    if Fields(3)==1; ATOM.altloc{end+1}=sscanf(A(i,17),'%c'); end	
    
    % Residue name.
    if Fields(4)==1; ATOM.resname{end+1}=sscanf(A(i,18:20),'%s'); end
    
    % Chain identifier.
    if Fields(5)==1; ATOM.chain{end+1}=sscanf(A(i,22),'%c'); end	
    
    % Residue sequence number.
    if Fields(6)==1; ATOM.resid(end+1)=sscanf(A(i,23:26),'%i'); end	
							  
    % Code for insertion of residues.      
    if Fields(7)==1; ATOM.icode{end+1}=sscanf(A(i,27),'%c'); end	
    
    % Orthogonal coordinates for X in Angstroms.
    if Fields(8)==1; ATOM.x(end+1)=sscanf(A(i,31:38),'%f'); end	
    
    % Orthogonal coordinates for Y in Angstroms.
    if Fields(9)==1; ATOM.y(end+1)=sscanf(A(i,39:46),'%f'); end	
    
    % Orthogonal coordinates for Z in Angstroms.
    if Fields(10)==1; ATOM.z(end+1)=sscanf(A(i,47:54),'%f'); end	
    
    % Occupancy
    if Fields(11)==1; ATOM.occupancy(end+1)=sscanf(A(i,55:60),'%f'); end
    
    % Temperature factor.    
    if Fields(12)==1; ATOM.beta(end+1)=sscanf(A(i,61:66),'%f'); end	
    
    % Segment identifier, left-justified.    
    if Fields(13)==1; ATOM.segname{end+1}=sscanf(A(i,73:76),'%s'); end
      
    %Element symbol, right-justified.      
    if Fields(14)==1; ATOM.element{end+1}=sscanf(A(i,77:78),'%s'); end

    % Charge on the atom.    
    if Fields(15)==1; ATOM.charge{end+1}=sscanf(A(i,79:80),'%s'); end
  end
  
  %Progress report block - optional
  CurrentProgressTime=toc;
  if (CurrentProgressTime-ProgressReportTime)>0.5	%Time to report progress
    try ProgressBarLocal(i,numlines);%Set current value on the local progress bar
    end
    ProgressReportTime=CurrentProgressTime;
  end	%Time to report progress
  
end

%disp([num2str(size(ATOM.index)), ' atoms structure read']);
%plot(ATOM.x)%!

%unzip_AGA(fname1,TempZipFileName);  %Remove Temporary Unzip file, if necesary

end

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

function [H,coords] = hessCalpha_mg(atoms,afftyPar,caidx) 
% set up the hessian matrix matrix between residues
% using ONLY the C-alpha atoms

  if nargin < 2
    afftyPar.minDist = 1.7; 
    afftyPar.maxDist = 12.0; 
    afftyPar.debug   = 1;
    afftyPar.dim     = 3;
  end
  minD2 = afftyPar.minDist^2;
  maxD2 = afftyPar.maxDist^2;
  % new code to get Hessian for 2d geometries
  dim   = afftyPar.dim;
  
  % which atoms are CA (C alphas)
%   caidx = find(strcmp(atoms.name,'CA'));
  
  % for ca residues that have alternate
  % coordinates deposited, pick only 
  % one set
  %[cai,caj] = unique(atoms.resid(caidx));
  %caidx = caidx(caj);
  
  % total # of C-alpha atoms 
  caLen = length(caidx);
  % where is the residue boundary
  aaidx = caidx - 1;

  % get all C-alpha atoms into one array
  if dim == 3
    coords = [atoms.x(caidx)' atoms.y(caidx)' atoms.z(caidx)'];
  else
    coords = [atoms.x(caidx)' atoms.y(caidx)'];
  end
  
  coords = coords';
  % how many of them in all?  
  atLen  = size(coords,2);
  % residue number look up table
  resId = 1:atLen;

  % hessian matrix
  H = sparse(dim*caLen,dim*caLen);
  %tmpHess = zeros(3,3);
  
  % distance computation 
  %tic;
%   size(coords,2)
%   coords

  for id = 1:size(coords,2)
      id;
    tmp = [];
    % coords is an array: dim X N
    % stand at an atom that is: coords(:,id)
    % and find its distance to all other atoms
    % look ahead in the computation, that is avoid
    % comparing with atoms that have already been looked at
    tmp = sum((coords(:,id+1:end) - repmat(coords(:,id),1,...
					 size(coords(:,id+1:end),2))).^2,1);
    % which of these distances satisfy the cutoff constraints
    idx = find((tmp >= minD2) & (tmp <= maxD2));
    %keyboard;
    % 
    tmp;
    idx;
    if idx
      % save the distances
      dist = tmp(idx);
      % since we look ahead, the find operator is unaware.
      % so add the location information.
      %idx = idx + id - 1;
      idx = idx + id ;
      %tmpHess(1:3,1:3) = 0;
      for i = 1:length(idx)
	c2 = coords(:,idx(i)) - coords(:,id);
	cc = -c2*c2'/dist(i);
	%tmpHess = tmpHess + cc;
	% super-block ij
	H((id-1)*dim+(1:dim),(idx(i)-1)*dim+(1:dim)) = cc;
	% super-block ji	
	H((idx(i)-1)*dim+(1:dim),(id-1)*dim+(1:dim)) = cc;	
	% diagonal-block ii
	H((id-1)*dim+(1:dim),(id-1)*dim+(1:dim)) = ...
	    H((id-1)*dim+(1:dim),(id-1)*dim+(1:dim)) -cc;      
	% diagonal-block jj
	H((idx(i)-1)*dim+(1:dim),(idx(i)-1)*dim+(1:dim)) = ...
	    H((idx(i)-1)*dim+(1:dim),(idx(i)-1)*dim+(1:dim))-cc;      
      end
      % superblock ii
      %H((id-1)*3+(1:3),(id-1)*3+(1:3)) = -tmpHess;      
    end
  end
end

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

function rc = close_dcd(h)

% rc = close_dcd(handle)

rc = fclose(h.fid);

end

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

function h = read_dcdheader(filename)

% new_handle = read_dcdheader(filename)

fid = fopen(filename, 'r', 'b');
i = fread(fid, 1, 'int32');
if i ~= 84              % wrong endianism; try again
  fclose(fid);
  fid = fopen(filename, 'r', 'l');
  i = fread(fid, 1, 'int32');
  if i ~= 84
    fclose(fid);
  end
end

h.fid = fid;

%
% Figure out how big the file is
%
fseek(fid,0,1);
h.endoffile = ftell(fid);
fseek(fid,4,-1);

% Read CORD 
s = fread(fid, 4, 'char');

% Read NSET
h.NSET = fread(fid, 1, 'int32');

% read ISTART
h.ISTART = fread(fid, 1, 'int32');

% read NSAVC
h.NSAVC = fread(fid, 1, 'int32');

% Reposition to 40 from beginning; read namnf, number of free atoms
fseek(fid, 40, -1);
h.NAMNF = fread(fid, 1, 'int32');

% Figure out if we're CHARMm or not
fseek(fid, 84, -1);
i = fread(fid, 1, 'int32');
if i == 0
  h.charmm = 0;
else
  h.charmm = 1;
  h.charmm_extrablock = 0;
  h.charmm_4dims = 0;
  % check for extra block
  fseek(fid, 48, -1);
  i = fread(fid, 1, 'int32');
  if i == 1
    h.charmm_extrablock = 1;
  end
  % check for 4 dims
  fseek(fid, 52, -1);
  i = fread(fid, 1, 'int32');
  if i == 1
    h.charmm_4dims = 1;
  end
end



% Read the timestep, DELTA.  Read a float if charmm, otherwise double
fseek(fid, 44, -1);
if h.charmm == 0
  h.DELTA = fread(fid, 1, 'float64');
else
  h.DELTA = fread(fid, 1, 'float32');
end

% Get the size of the next block, and skip it
% This is the title
fseek(fid, 92, -1);
newsize = fread(fid, 1, 'int32');
numlines = fread(fid, 1, 'int32');
fseek(fid, numlines*80, 0);
newsize = fread(fid, 1, 'int32');

% Read in a 4, then the number of atoms, then another 4
i = fread(fid, 1, 'int32');
h.N = fread(fid, 1, 'int32');
i = fread(fid, 1, 'int32');


% stuff with freeindexes.  Just smile and nod.
if h.NAMNF ~= 0
  fsize = fread(fid, 1, 'int32');  % should be N-NAMNF*4
  h.FREEINDEXES = fread(fid, h.N - h.NAMNF, 'int32');
  fsize = fread(fid, 1, 'int32');  % should be N-NAMNF*4
end 

end
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

function [x,y,z] = read_dcdstep(h)

%
% [x,y,z] = read_dcdstep(handle)
%

% If this is a CHARMm file and contains an extra data block, we must skip it
if h.charmm & h.charmm_extrablock
  blocksize = fread(h.fid, 1, 'int32');
  fseek(h.fid, blocksize, 0);
  blocksize = fread(h.fid, 1, 'int32');
end

if h.NAMNF == 0

  % Get x coordinates 
  blocksize = fread(h.fid, 1, 'int32');
  x = fread(h.fid, blocksize/4, 'float32');
  blocksize = fread(h.fid, 1, 'int32');

  % Get y coordinates 
  blocksize = fread(h.fid, 1, 'int32');
  y = fread(h.fid, blocksize/4, 'float32');
  blocksize = fread(h.fid, 1, 'int32');

  % Get z coordinates 
  blocksize = fread(h.fid, 1, 'int32');
  z = fread(h.fid, blocksize/4, 'float32');
  blocksize = fread(h.fid, 1, 'int32');

else    
 
  % this is not implemented in the VMD code I copied from  

end  
 
% Skip the 4th dimension, if present
if h.charmm & h.charmm_4dims
  blocksize = fread(h.fid, 1, 'int32');
  fseek(h.fid, blocksize, 0);
  blocksize = fread(h.fid, 1, 'int32');
end

end
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

function xyz = readdcd(filename, ind)

% xyz = readdcd(filename, indices)
% reads an dcd and puts the x,y,z coordinates corresponding to indices 
% in the rows of x,y,z

h = read_dcdheader(filename)
nsets = h.NSET;
natoms = h.N;
numind = length(ind);

x = zeros(natoms,1);
y = zeros(natoms,1);
z = zeros(natoms,1);

if nsets == 0
  xyz = zeros(1,3*numind);
  nsets = 99999;
else
  xyz = zeros(nsets, 3*numind);
end

for i=1:nsets
  pos = ftell(h.fid);
  if pos == h.endoffile 
    break;
  end
  [x,y,z] = read_dcdstep(h);
  xyz(i,1:3:3*numind) = x(ind)';
  xyz(i,2:3:3*numind) = y(ind)';
  xyz(i,3:3:3*numind) = z(ind)';
end

close_dcd(h);

end
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

function h = write_dcdheader(filename, nsets, natoms )

%
% handle = write_dcdheader(filename, nsets, natoms ) 
%

fid = fopen(filename, 'w');
h.fid = fid;

%
% write an 84 byte block header
%

fwrite(fid, 84, 'int32');             
fwrite(fid, 'CORD', 'uchar');         %  4
fwrite(fid, nsets, 'int32');          %  4
fwrite(fid, 0, 'int32');              %  4 - ISTART
fwrite(fid, 1, 'int32');              %  4 - NSAVC
fwrite(fid, zeros(6,1), 'int32');     % 24
fwrite(fid, 1, 'float64');            %  8 - DELTA
fwrite(fid, zeros(9,1), 'int32');   % 36
fwrite(fid, 84, 'int32');

%
% Now write a 164 byte block for the title
%

fwrite(fid, 164, 'int32'); 
fwrite(fid, 2, 'int32');     % number of 80-character lines in the title
count = fprintf(fid, 'REMARKS FILENAME=%s CREATED BY MATLAB',filename);
for i = count+1:80
  fwrite(fid, ' ', 'uchar');
end
count = fprintf(fid, 'REMARKS DATE: %s ', datestr(now));
for i = count+1:80
  fwrite(fid, ' ', 'uchar');
end
fwrite(fid, 164, 'int32');

%
% write the number of atoms
%

fwrite(fid, 4, 'int32');
fwrite(fid, natoms, 'int32');
fwrite(fid, 4, 'int32');

rc = 0;


end
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

function rc = write_dcdstep(h, x, y, z)

%
% rc = write_dcdstep, handle, x, y, z);
%

fwrite(h.fid, 4*length(x), 'int32');
fwrite(h.fid, x, 'float32');
fwrite(h.fid, 4*length(x), 'int32');

fwrite(h.fid, 4*length(x), 'int32');
fwrite(h.fid, y, 'float32');
fwrite(h.fid, 4*length(x), 'int32');

fwrite(h.fid, 4*length(x), 'int32');
fwrite(h.fid, z, 'float32');
fwrite(h.fid, 4*length(x), 'int32');

rc = 0;

end
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%

function rc = writedcd(filename, x, y, z)

%
% rc = writedcd(filename, x, y, z)
% x, y, z are of size [natoms nsets]
%

[natoms, nsets] = size(x);

h = write_dcdheader(filename, nsets, natoms);
for i=1:nsets
  rc = write_dcdstep(h, x(:,i), y(:,i), z(:,i));
end

close_dcd(h);
 
end

end

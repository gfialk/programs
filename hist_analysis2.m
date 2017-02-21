clear all;

% flag - chose if to save data of peaks as there values or in a binary
% form.
is_bin = 0;

base_dir = '../bam_to_window';
tissue_dir = dir([base_dir '/tissue*']);
num_of_tissues = length(tissue_dir);

WIN_SIZE = 10000;
CHROM_SIZE = 2;
WIND_CORD = 1;
NUM_READS = 3;
FILE_NAME = 1;



%initialize atlas metrix 

chrom_sizes = importdata('../programs/hg19_chrom_sizes.txt');
genome_size = int16(sum(chrom_sizes(:,CHROM_SIZE)/WIN_SIZE));
atlas_mat = zeros(num_of_tissues + 1, genome_size + 1);

index = 2;

for i=1:2 %%%%%%%%%%%%%%%%%%%  change to 1:25 %%%%%%%%%%%%%%%%%%%% 
    chr_size = chrom_sizes(i,CHROM_SIZE);
    for j=0:WIN_SIZE:chr_size
        % location is a catcatanation of chromosome number and cordination
        loc = str2num([num2str(i) num2str(j)]);
        atlas_mat(1, index) = loc;
        index = index+1;
    end
end


FILE_INDEX = 0;
%iterate through every modification directory
for tissue = tissue_dir'
    modification_files = dir([base_dir '/' tissue_dir.name '/*.win']);
    for modification = modification_files'
        
        FILE_INDEX = FILE_INDEX + 1;
        atlas_mat(FILE_INDEX + 1, FILE_NAME) = FILE_INDEX;
        fprintf('tissue: %s, modification: %s\n', tissue.name, modification.name);
        file = importdata([base_dir '/' tissue_dir.name '/' modification.name]);
        %fiealds of file are 'chromosome', 'start', 'end', 'number of reads'
        
        chr = file(:,1);
        data = file(:,2:4);
        
        % convert X,Y,M chromosomes to there int representation (23,24,25
        % respectively).    %I updated this into shell script and than can
        % erase these three lines
%         chr(:) = strrep(chr(:), 'X', '23');
%         chr(:) = strrep(chr(:), 'Y', '24');
%         chr(:) = strrep(chr(:), 'M', '25');
        
        % convert chromosome field (1st coulmn) to it's int representation
         %chr = cellfun(@str2num,chr(:));
        
               
        for k=1:2 %%%%%%%%%%%%%%%%%%%  change to 1:25 %%%%%%%%%%%%%%%%%%%%
            chr_k = find (chr == k);  
            chr_data = data(chr_k, :);
            
            % choose windows with number of read higher than the 95th
            % percentile
            cutoff = prctile(chr_data(:,NUM_READS), 90);
            
            
            for l=1:length(chr_k)
                loc = chr_data(l,WIND_CORD);
                cord = str2num([num2str(k) num2str(loc)]);
                num_reads = chr_data(l,NUM_READS);
                
                % binary option (represent windows with number of reads 
                % higher than cutoff - peaks - as 1)
                if num_reads >= cutoff & is_bin
                   atlas_mat(FILE_INDEX+1, atlas_mat(WIND_CORD,:) == cord) = 1;
                else
                    % no use of cutoff. report number of reads for every window 
                    atlas_mat(FILE_INDEX+1, atlas_mat(1,:) == cord) = num_reads;
                end
                                           
            end
        end
        
    end
end
                

%remove background signle by finding a cutoff of the 95'th percentile

%cutoff = prctile(num_reads, 95);


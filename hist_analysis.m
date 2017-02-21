clear all;
base_dir = 'bam_to_window';
win_files = dir([base_dir '/*.win']);
length(win_files)
for file = win_files'
    fprintf('file: %s\n', file.name);
    win_file = importdata(file.name);
    data = win_file.data;
    location = data(:,2);
    num_reads = data(:,3);
    plot(location, num_reads)
       
end





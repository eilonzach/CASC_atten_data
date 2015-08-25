dbdir = '/Users/zeilon/Work/CascAtten/CAdb';
dbnam = 'cascattendb';

run('/Users/zeilon/Documents/MATLAB/CascAtten/CascMAP.m')

ista = 1;
yrfiles = dir('dbdata/*.mat');
for iy = 1:length(yrfiles)
    load(sprintf('dbdata/%s',yrfiles(iy).name));
    yrdat = srStation;
    nyrstas = length(yrdat.name);
    
    

	figure(gcf), hold on
    scatter(yrdat.longitude,yrdat.latitude,65,'r','o','filled')
        
        
end
   return 

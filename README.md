# LiQCS
LiDAR Quality Control Suite
Created by Jordan Godau, Graeme Prendergast, Brett Edwards

LiQCS is intended to be an all-in-one lidar quality control suite.

Proper usage from the command line looks as follows:
    
    >> python liqcs.py -i C:\User\Lidar
        This command will run all QC checks on every .las/.laz file
        in the 'Lidar' folder. Output will be stored in the same folder.

    >> python liqcs.py -i C:\User\Lidar -o C:\User\Results
        This command will run all QC checks, then store output in
        the 'Results' folder.

    >> python liqcs.py -i C:\User\Lidar -t /i /n /dg /t
        This command runs only 4 QC checks, defined by each character
        in '/i /n /dg /t', and stores the results in the same folder.

*Each QC test has a shorthand to be defined after the -t flag as follows:
    /all    : run all tests
    /l      : lasinfo, generate lasinfo text files for each input
    /n      : check filename, verify if each filename adheres to the GeoBC specs
    /dl     : last density grid, create a last density grid .tif for each input
    /dg     : ground density grid, create a ground density grid .tif for each input
    /i      : intensity grid, create an intensity grid .tif for each input
    /t      : tile index, generate the bounding box tileset for the input batch
    /s      : lidary summary, summarize lasinfo outputs into a single csv file
    /h      : density histogram, plot point density histogram for the input batch
    /v      : void check, identify non-water point voids in the lidar, if any

*For example, consider the following command:
    *>> python liqcs.py -i C:\User\Lidar -t /dg /t /h
        This command would generate ground density grids, a tile index for the batch,
        and a density histogram for the batch.
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
**DOCUMENTATION WRITTEN BY GRAEME PRENDERGAST**

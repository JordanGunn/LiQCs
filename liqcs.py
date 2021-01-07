from lidar_summary import LidarSummary, find_format_errors
from gdalnumeric import BandWriteArray, BandReadAsArray, CopyDatasetInfo
from Create_Mask import create_mask
from functools import partial
from scipy.stats import norm
from laspy.file import File
from math import floor,ceil

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import multiprocessing as mp
import geopandas as gpd
import pandas as pd
import numpy as np

import concurrent.futures
import subprocess
import argparse
import sys
import os

import gdal
import osr
import shapely

import struct
import glob
import csv
import re

import datetime
import time

def main():
    """
    Parse and interpret command line arguments, and then run all tests specified
    by the user. Manage control flow between I/O and test functions.
    """
    arg_parser = argparse.ArgumentParser(
        prog='liqcs.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=
r"""
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

Each QC test has a shorthand to be defined after the -t flag as follows:
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

For example, consider the following command:
    >> python liqcs.py -i C:\User\Lidar -t /dg /t /h
        This command would generate ground density grids, a tile index for the batch,
        and a density histogram for the batch.
"""
)
    arg_parser.add_argument('-i', '--input', help="input directory path (ie, C:\\User\\...)")
    arg_parser.add_argument('-fl', '--filelist', help="text file containing a list of files seperated by newlines",
                                default=None, nargs=1, type=str, metavar="<<file_list.txt>>")
    arg_parser.add_argument('-o', '--output', help="output directory path (default is input dir)")
    arg_parser.add_argument('-t', '--tests', help="which QC tests to run (default is /h /n /t /s /l)",
                                default=r'/h /n /t /s /l', nargs='+', type=str)
    arg_parser.add_argument('-c', '--cores', help="how many processor cores to use",
                                type=int, default=4)
    arg_parser.add_argument('-s', '--shapefile', help='shapefile paths for AOI and breaklines. Needed if using /v (void detection) test flag',
                                 nargs=2, type=str, default=None, metavar=('<breaklines.shp>', '<AOI.shp>'))
    args = arg_parser.parse_args()

    # Defining an input directory path
    if args.input == None and args.filelist == None:
        sys.exit("You must specify an input directory! Use [-h] for help")

    if args.input is not None and args.filelist is None:   
        input_path = args.input

        # Defining an output directory path
        if args.output == None:
            results_path = f'{input_path}\\liqcs_results'
        else:
            results_path = args.output
        if not os.path.exists(results_path):
            os.makedirs(results_path)

        # Using glob to gather all the .las or .laz files in the input directory
        infile_glob = [f for f in glob.glob(f'{input_path}\\**\*',recursive=True)
                    if f[-4:] in {'.las','.laz'}]
        
        n = len(infile_glob)
        if n == 0:
            sys.exit("Couldn't find any .las or .laz files in the input directory")

    if args.input is None and args.filelist is not None:   
        sep = '\\'
        input_path = args.filelist[0].split(sep)
        input_path.remove(input_path[-1])
        input_path = sep.join(input_path)
        print(args.filelist[0])

        # Defining an output directory path
        if args.output == None:
            results_path = f'{input_path}\\liqcs_results'
        else:
            results_path = args.output
        if not os.path.exists(results_path):
            os.makedirs(results_path)

        with open(args.filelist[0], 'r') as f:
            infile_glob = [line.strip() for line in f]
        n = len(infile_glob)
        if n == 0:
            sys.exit("Couldn't find any .las or .laz files in the input directory")          


    # Set flags for which tests are to be run
    test_flags, epsg_code = set_test_flags(args.tests)

    # Define a shapefile; we're totally fine with this being None
    shapefile = args.shapefile

    start_time = time.time()

    if '/n' in test_flags: # Check for bad filenames
        # Create a subdirectory for format reports
        format_path = f'{results_path}\\formatting'
        if not os.path.exists(format_path):
            os.makedirs(format_path)
        print(f"Checking {n} filenames...")
        bad_filenames = [f.split('\\')[-1] for f in infile_glob
                         if filename_is_bad(f)]

        with open(f'{format_path}\\filename_check.txt', mode='w') as out:
            print(f"Bad filenames: {len(bad_filenames)} of {len(infile_glob)}", file=out)
            print(f"\nList of bad filenames:", file=out)
            for fn in bad_filenames:
                print(fn, file=out)

    if '/t' in test_flags:
        # Running LiDAR_TileIndex_GPKG (kinda) to generate a tileset for the lidar
        print("Generating tile index...")
        tileset = generate_tile_index(infile_glob,epsg_code)
        tileset.to_file(f'{results_path}\\bbox_tiles_delivered_{datetime.date.today()}.gpkg', driver='GPKG')


    if '/l' in test_flags: # Generate lasinfo.txt's
        # Create a subdirectory for the lasinfo.txts
        lasinfo_path = f'{results_path}\\lasinfo'
        if not os.path.exists(lasinfo_path):
            os.makedirs(lasinfo_path)
            os.makedirs(f'{lasinfo_path}\\last_return')
            os.makedirs(f'{lasinfo_path}\\ground')

        # Make a list of files .txt to manage extra large input directories
        lof_path = f'{input_path}\\liqcs_lof.txt'

        # Loop through infile glob, generating two lasinfo.txts for each
        print(f"Generating {n} lasinfos...")
        os.environ["PATH"] += os.pathsep + "C:\\LAStools\\bin"

        with open(lof_path, mode='w') as out:
            for fp in infile_glob:
                fn = fp.split('\\')[-1].split('.')[0]
                if not os.path.exists(f"{lasinfo_path}\\last_return\\{fn}.txt"):
                    print(f'{fp}', file=out)

        subprocess.call(
        f'lasinfo.exe -lof "{lof_path}" -cd -otxt -odir {lasinfo_path}\\last_return '
        f'-cores {args.cores}', shell=True
        )
        os.remove(lof_path)

        with open(lof_path, mode='w') as out:
            for fp in infile_glob:
                fn = fp.split('\\')[-1].split('.')[0]
                if not os.path.exists(f"{lasinfo_path}\\ground\\{fn}_ground.txt"):
                    print(f'{fp}', file=out)

        subprocess.call(
        f'lasinfo.exe -lof "{lof_path}" -cd -otxt -odix _ground -keep_class 2 '
        f'-odir {lasinfo_path}\\ground -cores {args.cores}', shell=True
        )
        os.remove(lof_path) # Get rid of the lof.txt once we're done with it


    if '/s' in test_flags or '/h' in test_flags:
        # Create a subdirectory for format reports
        format_path = f'{results_path}\\formatting'
        if not os.path.exists(format_path):
            os.makedirs(format_path)
        # Obviously, this should only be done if lasinfos have been generated
        lasinfo_lr_glob = glob.glob(f'{lasinfo_path}\\last_return\\*.txt')
        lasinfo_g_glob = glob.glob(f'{lasinfo_path}\\ground\\*.txt')

        if '/s' in test_flags:
            # Running lidar_summary() to generate a summary .csv from lasinfo.txts
            print("Generating lidar summary...")
            for file in lasinfo_lr_glob:
                info_record = LidarSummary(file)
                summary_csv = info_record.output(format_path, write_errors=False)

            find_format_errors(summary_csv, format_path, save_output=True)

        if '/h' in test_flags:
            # Running calculate_density to generate a density histograms
            print("Generating density histograms...")
            make_density_histograms(lasinfo_lr_glob, lasinfo_g_glob, results_path)

    if '/dg' in test_flags or '/dl' in test_flags or '/i' in test_flags:
        # Create subdirectories to store the grids
        grid_path = f'{results_path}\\grids'
        if not os.path.exists(grid_path):
            os.makedirs(grid_path)
        if '/i' in test_flags:
            if not os.path.exists(f'{grid_path}\\intensity'):
                os.makedirs(f'{grid_path}\\intensity')
        if '/dg' in test_flags or '/dl' in test_flags:
            if not os.path.exists(f'{grid_path}\\density'):
                os.makedirs(f'{grid_path}\\density')
            if '/dg' in test_flags:
                if not os.path.exists(f'{grid_path}\\density\\ground_density'):
                    os.mkdir(f'{grid_path}\\density\\ground_density')
            if '/dl' in test_flags:
                if not os.path.exists(f'{grid_path}\\density\\last_density'):
                    os.mkdir(f'{grid_path}\\density\\last_density')

        print(f"Generating {n} sets of grids...")

        # We want to prioritize doing 2 at a time to minimize decompression time
        if '/dg' in test_flags and '/dl' in test_flags and '/i' in test_flags:
            generate_grids_part = partial(generate_grids, type='all', grid_path=grid_path, epsg_code=epsg_code)  
            with concurrent.futures.ProcessPoolExecutor(args.cores) as executor:
                # args_pack = [('intensity',fp,grid_path,epsg_code) for fp in infile_glob]
                # pool.starmap(generate_grids, args_pack)
                executor.map(generate_grids_part, infile_glob)

        # if ground and last
        elif '/dg' in test_flags and '/dl' in test_flags and '/i' not in test_flags:
            generate_grids_part = partial(generate_grids, type=['ground', 'last'], grid_path=grid_path, epsg_code=epsg_code)  
            with concurrent.futures.ProcessPoolExecutor(args.cores) as executor:
                # args_pack = [('intensity',fp,grid_path,epsg_code) for fp in infile_glob]
                # pool.starmap(generate_grids, args_pack)
                executor.map(generate_grids_part, infile_glob)
        
        # if intensity and last
        elif '/i' in test_flags and '/dl' in test_flags and '/dg' not in test_flags:
            generate_grids_part = partial(generate_grids, type=['intensity', 'last'], grid_path=grid_path, epsg_code=epsg_code)  
            with concurrent.futures.ProcessPoolExecutor(args.cores) as executor:
                # args_pack = [('intensity',fp,grid_path,epsg_code) for fp in infile_glob]
                # pool.starmap(generate_grids, args_pack)
                executor.map(generate_grids_part, infile_glob)

        # if intensity and ground
        elif '/i' in test_flags and '/dg' in test_flags and '/dl' not in test_flags:
            generate_grids_part = partial(generate_grids, type=['intensity', 'ground'], grid_path=grid_path, epsg_code=epsg_code)  
            with concurrent.futures.ProcessPoolExecutor(args.cores) as executor:
                # args_pack = [('intensity',fp,grid_path,epsg_code) for fp in infile_glob]
                # pool.starmap(generate_grids, args_pack)
                executor.map(generate_grids_part, infile_glob)

        # if only intensity grids
        elif '/i'  in test_flags:
            generate_grids_part = partial(generate_grids, type='intensity', grid_path=grid_path, epsg_code=epsg_code)  
            with concurrent.futures.ProcessPoolExecutor(args.cores) as executor:
                # args_pack = [('intensity',fp,grid_path,epsg_code) for fp in infile_glob]
                # pool.starmap(generate_grids, args_pack)
                executor.map(generate_grids_part, infile_glob)
        
        # if only ground density
        elif '/dg' in test_flags: 
            generate_grids_part = partial(generate_grids, type='ground', grid_path=grid_path, epsg_code=epsg_code)  
            with concurrent.futures.ProcessPoolExecutor(args.cores) as executor:
                # args_pack = [('intensity',fp,grid_path,epsg_code) for fp in infile_glob]
                # pool.starmap(generate_grids, args_pack)
                executor.map(generate_grids_part, infile_glob)

        # if only last density
        elif '/dl' in test_flags:
            generate_grids_part = partial(generate_grids, type='last', grid_path=grid_path, epsg_code=epsg_code)  
            with concurrent.futures.ProcessPoolExecutor(args.cores) as executor:
                # args_pack = [('intensity',fp,grid_path,epsg_code) for fp in infile_glob]
                # pool.starmap(generate_grids, args_pack)
                executor.map(generate_grids_part, infile_glob)
       

    if '/v' in test_flags:
        # Create a density grid glob
        density_glob = glob.glob(f'{grid_path}\\density\\last_density\\*.tif')
        # Create a subdirectory for the void grids
        void_path = f'{results_path}\\voids'
        if not os.path.exists(void_path):
            os.makedirs(void_path)
        # create mask for void detection routine
        create_mask(breaklines=shapefile[0], project_area=shapefile[1], outdir=void_path)
        # Loop through the infiles, generating void grids and summarizing results
        with open(os.path.join(results_path, 'voids_summary.csv'), 'a', newline='') as csv_file:
            # create function reference for loop optimization
            writer = csv.writer(csv_file, delimiter=',')
            # write header names to csv file
            writer.writerow(['filename', 'size of void (pixels)'])
            for i,density_grid in enumerate(density_glob):
                print(f"Generating void grids: {i+1} of {n}\r", end='')
                void = void_grid(density_grid, void_path, os.path.join(void_path, 'mask.shp')) # create void grid
                if void > 500:
                    row = os.path.basename(density_grid).replace('void_', ''), void
                    writer.writerow(row)
            print()


    print("\nRuntime: %f seconds" % (time.time()-start_time))


def set_test_flags(arg_string):
    """Return a complete list set of tests to be run, including dependencies"""
    # We use a set to check membership in O(1) time expected
    tests = set(arg_string)
    epsg = None

    # if all shortcut specified, add all tests
    if '/all' in tests:
        tests.add('/s')
        tests.add('/h')
        tests.add('/l')
        tests.add('/i')
        tests.add('/n')
        tests.add('/dg')
        tests.add('/dl')
        tests.add('/v')
        tests.add('/t')

    # else handle tests that depend on other tests      
    else:
        # Lidar Summary (s) and Density Histogram (d) have dependencies on lasinfo (i)
        if ('/s' in tests or '/h' in tests) and ('/l' not in tests):
            print("Specified a test that requires lasinfo results; "
                "adding lasinfo to tests")
            tests.add('/l')

        # Void check (v) requires a density grid (g)
        if ('/v' in tests and '/dl' not in tests):
            print("Specified a test that requires a density grid; "
                "adding density grid to tests")
            tests.add('/dl')

    # Other dependency logic goes here

    # Prompt user for an EPSG code if needed
    if ('/dl' in tests or '/dg' in tests or '/i' in tests or '/t' in tests):
        print(
            "Specified a test that requires defining a CRS\n"
            "\nEPSG 2955: NAD83(CSRS) / UTM zone 11"
            "\nEPSG 3157: NAD83(CSRS) / UTM zone 10"
            "\nEPSG 3156: NAD83(CSRS) / UTM zone 9"
            "\nEPSG 3155: NAD83(CSRS) / UTM zone 8"
            "\nEPSG 3154: NAD83(CSRS) / UTM zone 7"
            "\nEPSG 3005: NAD83 / BC Albers\n")
        epsg = int(input("Please enter EPSG code:\t"))
        print()

    return tests, epsg


def generate_grids(infile_path, type, grid_path, epsg_code):
    """Grid a given lidar file by density and/or intensity onto a .tif."""
    filename = infile_path.split('\\')[-1] # Define the infile name
    results_path = '\\'.join(grid_path.split('\\')[:-1])

    try:
        if type == 'all' and not (
                os.path.exists(f'{grid_path}\\density\\ground_density\\{filename[:-4]}_ground_density.tif')
                and os.path.exists(f'{grid_path}\\intensity\\{filename[:-4]}_intensity_grid.tif')
                and os.path.exists(f'{grid_path}\\density\\last_density\\{filename[:-4]}_last_density.tif')
            ): # Decompress once, generate all grids
            lidar_grid = Grid()
            lidar_grid.proj_from_epsg(epsg_code)
            lidar_grid.read_lidar(infile_path,results_path)
            
            lidar_grid.create_grid(last_grid=True, filter_by_class=2,class_grid=True, intensity_grid=True)
            lidar_grid.save_grid(gridtype=f'{filename[:-4]}_ground',
                                    outdir=f'{grid_path}\\density\\ground_density')
            lidar_grid.save_grid(gridtype=f'{filename[:-4]}_last',
                                    outdir=f'{grid_path}\\density\\last_density')
            lidar_grid.save_grid(gridtype=f'{filename[:-4]}_intensity',
                                    outdir=f'{grid_path}\\intensity')


        elif 'last' in type and 'ground' in type and not os.path.exists(f'{grid_path}\\density\\last_density\\{filename[:-4]}_last_grid.tif') and not os.path.exists(f'{grid_path}\\density\\ground_density\\{filename[:-4]}_ground_grid.tif'
            ): # Decompress once, generate a density grid
            lidar_grid = Grid()
            lidar_grid.proj_from_epsg(epsg_code)
            lidar_grid.read_lidar(infile_path,results_path)

            lidar_grid.create_grid(last_grid=True, filter_by_class=2, class_grid=True)
            lidar_grid.save_grid(gridtype=f'{filename[:-4]}_ground',
                                    outdir=f'{grid_path}\\density\\ground_density')
            lidar_grid.save_grid(gridtype=f'{filename[:-4]}_last',
                                    outdir=f'{grid_path}\\density\\last_density')

        elif 'last' in type and 'intensity' in type and not os.path.exists(f'{grid_path}\\density\\last_density\\{filename[:-4]}_last_grid.tif') and not os.path.exists(f'{grid_path}\\intensity\\{filename[:-4]}_intensity_grid.tif'
            ): # Decompress once, generate a density grid
            lidar_grid = Grid()
            lidar_grid.proj_from_epsg(epsg_code)
            lidar_grid.read_lidar(infile_path,results_path)

            lidar_grid.create_grid(last_grid=True, intensity_grid=True)
            lidar_grid.save_grid(gridtype=f'{filename[:-4]}_last',
                                    outdir=f'{grid_path}\\density\\last_density')
            lidar_grid.save_grid(gridtype=f'{filename[:-4]}_intensity',
                                    outdir=f'{grid_path}\\intensity')

        elif 'ground' in type and 'intensity' in type and not os.path.exists(f'{grid_path}\\intensity\\{filename[:-4]}_intensity_grid.tif') and not os.path.exists(f'{grid_path}\\density\\ground_density\\{filename[:-4]}_ground_grid.tif'
            ): # Decompress once, generate a density grid
            lidar_grid = Grid()
            lidar_grid.proj_from_epsg(epsg_code)
            lidar_grid.read_lidar(infile_path,results_path)
            lidar_grid.create_grid(filter_by_class=2, class_grid=True, intensity_grid=True)
            lidar_grid.save_grid(gridtype=f'{filename[:-4]}_ground',
                                    outdir=f'{grid_path}\\density\\ground_density')
            lidar_grid.save_grid(gridtype=f'{filename[:-4]}_intensity',
                                    outdir=f'{grid_path}\\intensity')

        elif type == 'last' and not os.path.exists(
            f'{grid_path}\\density\\last_density\\{filename[:-4]}_last_grid.tif'
            ): # Decompress once, generate a density grid
            lidar_grid = Grid()
            lidar_grid.proj_from_epsg(epsg_code)
            lidar_grid.read_lidar(infile_path,results_path)

            lidar_grid.create_grid(last_grid=True)
            lidar_grid.save_grid(gridtype=f'{filename[:-4]}_last',
                                    outdir=f'{grid_path}\\density\\last_density')

        elif type == 'ground' and not os.path.exists(
                f'{grid_path}\\density\\ground_density\\{filename[:-4]}_ground_grid.tif'
            ): # Decompress once, generate a density grid
            lidar_grid = Grid()
            lidar_grid.proj_from_epsg(epsg_code)
            lidar_grid.read_lidar(infile_path,results_path)

            lidar_grid.create_grid(filter_by_class=2, class_grid=True)
            lidar_grid.save_grid(gridtype=f'{filename[:-4]}_ground',
                                    outdir=f'{grid_path}\\density\\ground_density')     

        elif type == 'intensity' and not os.path.exists(
                f'{grid_path}\\intensity\\{filename[:-4]}_intensity_grid.tif'
            ): # Decompress once, generate an intensity grid
            lidar_grid = Grid()
            lidar_grid.proj_from_epsg(epsg_code)
            lidar_grid.read_lidar(infile_path,results_path)

            lidar_grid.create_grid(intensity_grid=True)
            lidar_grid.save_grid(gridtype=f'{filename[:-4]}_intensity',
                                    outdir=f'{grid_path}\\intensity')

    except Exception as e:
        print(f'{filename} could not be written to a grid')
        with open(f'{results_path}\\grid_failures.txt', 'a+') as f:
            f.write(f'{filename} :: {e}\n')
        pass


def make_density_histograms(lasinfo_lr_glob, lasinfo_g_glob, outdir):
    """Creating last return and ground point density histograms"""
    try:
        lr_dens_none = [lasinfo_density(f, last_return=True) for f in lasinfo_lr_glob]
        lr_dens = [d for d in lr_dens_none if d]  # remove values equal to None   !!! Fixed bug causing ground dens histogram to not be created 2019-10-24 :: By Jordan !!!
        density_histogram(lr_dens,outdir+'\\last_return_density')
    except:
        print('Error creating last return density statistics')
        pass

    try:
        ground_dens_none = [lasinfo_density(f, last_return=False) for f in lasinfo_g_glob]  # may have values equal to None
        ground_dens = [d for d in ground_dens_none if d]  # remove values equal to None 
        density_histogram(ground_dens,outdir+'\\ground_density')
    except:
        print('Error creating ground density statistics, file may contain no ground points')
        pass


def void_grid(file, outdir, mask=None):
    """Detect voids in a density grid .tif, and output a void grid and csv row."""
    fn = file.split('\\')[-1][:-14]
    out_path = f'{outdir}\\{fn}_voids.tif'

    grid1 = gdal.Open(file, gdal.GA_ReadOnly) # Read the density grid

    band = grid1.GetRasterBand(1) # extract raster band from density grid
    wkt = grid1.GetProjection() # extract projection from density grid
    band_array = BandReadAsArray(band) # read raster band as an array

    # Anywhere density > 0 is set to 1, else unchanged
    dataOut = np.where(band_array > 0, 1, band_array)
    dataOut[0:20] = 1
    dataOut[-20::] = 1
    dataOut[:, 0:20] = 1
    dataOut[:, -20::] = 1

    driver = gdal.GetDriverByName("GTiff") # Create Geotiff driver object

    # Create new geotiff file object with same dimensions and data-type as input file
    dsOut = driver.Create(out_path, grid1.RasterXSize, grid1.RasterYSize, 1,band.DataType)

    CopyDatasetInfo(grid1, dsOut) # copy properties of density grid to new output raster

    dsOut.SetProjection(wkt) # set projection for new raster
    # assign band attribute to varaible
    bandOut = dsOut.GetRasterBand(1)  # set number of bands for new raster

    # bandOut.SetNoDataValue(0)
    # write new information to blank raster
    BandWriteArray(bandOut, dataOut)
    gdal.SieveFilter(bandOut, None, bandOut, 25, connectedness=8)
    # rasterize polygon to array, burn values
    if mask is not None:
        gdal.Rasterize(dsOut, mask, burnValues=[1])
    
    # filter out peppery void areas smaller than 9x9 pixels 
    # gdal.SieveFilter(bandOut, None, bandOut, 16, connectedness=8)
    sieve_band = dsOut.GetRasterBand(1)
    sieve_array = BandReadAsArray(sieve_band)
    
    # calcuate number of pixels missing from the raster
    void_area = int(np.size(sieve_array) - np.sum(sieve_array))

    # flush data to output
    bandOut.FlushCache()

    # close the datasource
    grid1 = None
    dsOut = None

    return void_area


def lasinfo_density(lasinfo, last_return=False):
    """Extract density from lasinfo txt file"""
    # open lasinfo file
    open_info = open(lasinfo, "r")
    # find index value of point density information
    for line in open_info:
        if 'point density' in line:
            # get density values and convert to float
            try:
                if last_return:
                    dens = float(line.split(' ')[-4])
                else:
                    dens = float(line.split(' ')[4])
                return dens
            except:
                print(f'Error getting density value from {lasinfo}')
                pass
    
    open_info.close()


def density_histogram(dataset, path, custom_title=None):
    """create histogram for density values"""
    # best fit of data
    (mu, sigma) = norm.fit(dataset)

    # the histogram of the data
    n, bins, patches = plt.hist(dataset, 15, density=1, facecolor='blue', alpha=0.75)

    # add a 'best fit' line
    # Changed mlab.normpdf() to scipy.stats.norm.pdf(), mlab version is deprecated
    y = norm.pdf(bins, mu, sigma)
    l = plt.plot(bins, y, 'y--', linewidth=2)

    #plot
    plt.xlabel('Density (points/sq. m)')
    plt.ylabel('Frequency')
    if custom_title is not None:
            plt.title = (f'{custom_title}:  mu={mu}  sigma={sigma}')
    else:
        plt.title(r'$\mathrm{Histogram\ of\ Density:}\ \mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))

    plt.grid(True)

    plt.savefig(f'{path}.png')
    plt.close()


def filename_is_bad(fp):
    """Check if a lidar filename is "bad," aka does not fit GeoBC naming conventions."""
    fn_regex = (
    'bc_[0-9]{3}[a-z][0-9]{3}_([0-9]_){3}x(yes|no)_'
    '[0-9]{1,2}_(bcalb|(utm[0-9]{1,2}))_[0-9]+\.la(s|z)'
    )
    fn = fp.split('\\')[-1]
    return not bool(re.fullmatch(fn_regex,fn))


def generate_tile_index(infile_list,epsg):
    """
    Generate a GeoPackage containing a set of tiles, where each tile is a
    2d bounding box around one of the input lidar point clouds
    """
    gdf = create_polygon_gdf(epsg_code=epsg)
    tiles = [create_bounding_box(f) for f in infile_list]

    return append_to_gdf(tiles,infile_list,gdf)


def create_bounding_box(file):
    """Create a bounding box from LiDAR data and return WKT polygon"""
    # parse las/laz header body
    hdr = parse_header(file, verbose=False)

    # gt min/max from las/laz header
    xmin, xmax, ymin, ymax = [hdr['xmin'], hdr['xmax'], hdr['ymin'], hdr['ymax']]

    # Create polygon from bounding box of las file
    bbox_poly = shapely.geometry.Polygon(
        [
            (xmin, ymax),
            (xmax, ymax),
            (xmax, ymin),
            (xmin, ymin),
            (xmin, ymax),
        ]
    )

    return bbox_poly


def create_polygon_gdf(epsg_code=None):
    """Create a GeoDataFrame for a polygon, and define CRS if possible."""
    # create empty geopandas dataframe and assign projection
    df = gpd.GeoDataFrame(columns=["map_tile", "file", "geometry"])

    # assign spatial reference system metadata
    if epsg_code is not None:
        df.crs = {"init": f"epsg:{epsg_code}"}

    return df


def append_to_gdf(geom_list, file_list, gdf):
    """Append a list of geometry in WKT format to a GeoDataFrame"""
    for geom, file in zip(geom_list, file_list):

        filename = os.path.basename(file)
        tile = filename[3:16]

        # append to dataframe
        gdf = gdf.append(
            {"map_tile": tile, "file": file, "geometry": geom},
            ignore_index=True
        )

    return gdf


def parse_header(filename, verbose=False):
    """Parse a las/laz file's header into a workable struct format."""
    headerstruct = (
        ('filesig', 4,'c',4) ,
        ('filesourceid' , 2,'H',1) ,
        ('reserved'     , 2,'H',1) ,
        ('guid1'        , 4,'L',1) ,
        ('guid2'        , 2,'H',1) ,
        ('guid3'        , 2,'H',1) ,
        ('guid4'        , 8,'B',8) ,
        ('vermajor'     , 1,'B',1) ,
        ('verminor'     , 1,'B',1) ,
        ('sysid'        , 32,'c',32) ,
        ('gensoftware'  , 32,'c',32) ,
        ('fileday'      , 2,'H',1) ,
        ('fileyear'     , 2,'H',1) ,
        ('headersize'   , 2,'H',1) ,
        ('offset'       , 4,'L',1) ,
        ('numvlrecords' , 4,'L',1) ,
        ('pointformat'  , 1,'B',1) ,
        ('pointreclen'  , 2,'H',1) ,
        ('numptrecords' , 4,'L',1) ,
        ('numptbyreturn', 20,'L',5) ,
        ('xscale'       , 8,'d',1) ,
        ('yscale'       , 8,'d',1) ,
        ('zscale'       , 8,'d',1) ,
        ('xoffset'      , 8,'d',1) ,
        ('yoffset'      , 8,'d',1) ,
        ('zoffset'      , 8,'d',1) ,
        ('xmax'         , 8,'d',1) ,
        ('xmin'         , 8,'d',1) ,
        ('ymax'         , 8,'d',1) ,
        ('ymin'         , 8,'d',1) ,
        ('zmax'         , 8,'d',1) ,
        ('zmin'         , 8,'d',1) ,
        ('waveform'     , 8,'Q',1) ,
        ('firstEVLR'    , 8,'Q',1) ,
        ('numEVLR'      , 4,'L',1) ,
        ('exnumbptrec'  , 8,'Q',1) ,
        ('exnumbyreturn',120,'Q',15)
        )

    header = {'infile':filename}

    with open(filename, 'rb') as fh:
        for i in headerstruct:
            if i[2] == 'c':
                value = fh.read(i[1])
            elif i[3] > 1:
                value = struct.unpack( '=' + str(i[3]) + i[2] , fh.read(i[1]) )
            else:
                value = struct.unpack( '=' + i[2] , fh.read(i[1]) )[0]
            if verbose:
                print(i[0] + '\t', i[2] + '\t', value)

            header[i[0]] = value

    # if laz file, subtract 128 from the point data format
    # (laz compression adds 128 to the point format)
    if header['infile'].endswith('.laz'):
        header['pointformat'] = header['pointformat']-128

    return header


class Grid:
    """Class to generate density and/or intensity grids from lidar."""
    def __init__(self):
       self.path = None #
       self.file_object = None #
       self.cell_size = 1 #
       self.row = None #
       self.col = None #
       self.NoDATA = None
       self.crs = None #
       self.min = None
       self.max = None
       self.verbose = False


    def read_lidar(self, filepath, outdir):
        """Method for reading las or laz data with error handling."""
        # handle laz file if it is encountered
        # Open las/laz file in read mode

        self.path = filepath

        if 'laz' in os.path.splitext(filepath)[1]:
            lastools_path = os.pathsep + 'C:\\LAStools\\bin'
            os.environ["PATH"] += lastools_path
            if self.verbose:
                print('\nDecompressing LAZ file...\n')
            try:
                self.file_object = File(filepath, mode='r')
                if self.verbose:
                    print('Finished decompressing')

            except Exception as e:
                print(f'Error opening {filepath}, file may be corrupt...\n')
                with open(f"{outdir}\\lidar_read_errors.txt", mode='a+') as f:
                    f.write(f"{filepath}: {e}\n")
            env_var = os.environ["PATH"]
            env_var.replace(lastools_path, '')  
        else:
            try:
                self.file_object = File(filepath, mode='r')
            except Exception as e:
                print(f'Error opening {filepath}, file may be corrupt...\n')
                with open(f"{outdir}\\lidar_read_errors.txt", mode='a+') as f:
                    f.write(f"{filepath}: {e}\n")
        if self.verbose:
            print(f'Created laspy file object from {filepath}')
        return self.file_object


    def proj_from_epsg(self, EPSG_code):
        """Method to generate OGC wkt georeferencing metadata from EPSG code."""

        self.crs = None
        sr = osr.SpatialReference() # Setup georeferencing info
        sr.ImportFromEPSG(EPSG_code)  # pass in the epsg code
        self.crs = sr.ExportToWkt()  # export georeferencing metadata in Well-Known Text format
        return self.crs


    def create_density_grid(self, last_return=False, filter_by_class=None):
        """
        Method to create density grid from lidar. User must have grid.file_object attribute
        set first using the read_lidar method. User can filter by last return or class code.

        filter by keyword argument parameters, either filter by class, last return, both, or neither
        then extract x, y, z values from file to be used later for density calculations
        """
        self.NoDATA = 0  # !!!!!!!!!!!!!!! Changed from -1 to 0 on 2019-08-13 !!!!!!!!!!!!!!!!
        self.array = None
        self.grid_type = 'density'

        las = self.file_object
        # extract minimum and maximum from input file; floor min & ceiling max

        self.min = [floor(las.header.min[0]), floor(las.header.min[1])]
        # self.max = [ceil(las.header.max[0]), ceil(las.header.max[1])]
        self.max = [floor(las.header.max[0])+1, floor(las.header.max[1])+1]

        # get x and axis distance (m) from las file
        dist_x = self.max[0] - self.min[0]
        dist_y = self.max[1] - self.min[1]
        # calculate number of columns for raster grid
        self.col = int(dist_x/self.cell_size)
        # self.col += 1
        self.row = int(dist_y/self.cell_size)
        # self.row += 1

        if last_return and filter_by_class is None:
            las_return_x, las_return_y = (
                las.x[las.num_returns == las.return_num],
                las.y[las.num_returns == las.return_num]
                )
            return_filter = True
            class_filter = False
            filter_filter = False

        elif filter_by_class is not None and last_return is False:
            las_class_x, las_class_y = (
                las.x[las.Classification == filter_by_class],
                las.y[las.Classification == filter_by_class]
                )
            class_filter = True
            return_filter = False
            filter_filter = False

        elif filter_by_class is not None and last_return:
            las_filter_x, las_filter_y = (
                las.x[las.Classification == filter_by_class
                    and las.num_returns == las.return_num],
                las.y[las.Classification == filter_by_class
                    and las.num_returns == las.return_num]
                )
            filter_filter = True
            class_filter = False
            return_filter = False

        else:
            las_x, las_y = las.x, las.y
            class_filter = False
            return_filter = False
            filter_filter = False

        # Create empty numpy array to write values to
        count = np.zeros((self.row, self.col)).astype(np.int32)

        # Apply -1 to have negative y resolution for raster
        ycell = -1 * self.cell_size

        # Scale or "project" values  of LiDAR data to grid,
        if return_filter:
            scale_x = (las_return_x - self.min[0]) / self.cell_size
            scale_y = (las_return_y - self.min[1]) / ycell

        elif class_filter:
            scale_x = (las_class_x - self.min[0]) / self.cell_size
            scale_y = (las_class_y - self.min[1]) / ycell

        elif filter_filter:
            scale_x = (las_filter_x - self.min[0]) / self.cell_size
            scale_y = (las_filter_y - self.min[1]) / ycell
        else:
            scale_x = (las_x - self.min[0]) / self.cell_size
            scale_y = (las_y - self.min[1]) / ycell

        # change type to integer and save as variables to use for index values
        ind_x = scale_x.astype(np.int32)
        ind_y = scale_y.astype(np.int32)

        # Loop through LiDAR point records, count, and add to raster grid
        if self.verbose:
            print(f'\nCalculating density for {os.path.splitext(os.path.basename(self.path))[0]}...')

        # Runtime bottleneck - This is O(n) WRT the number of points in the point cloud
        # We can't do better than this withouth sacrificing accuracy
        for x, y in np.nditer([ind_x, ind_y]):
            count[y, x] += 1

        # Fill areas areas lacking data with keyword argument specified no data value
        count_noData = (np.where(count > 0, count, self.NoDATA)).astype(np.int32)

        # calculate density
        self.array = (count_noData / self.cell_size).astype(np.int32)
        # lazy fix for top row of cells being generated outside AOI inexplicably
        self.array = np.delete(self.array, 0, axis=0)


    def create_intensity_grid(self):
        """
        Method to create intensity grid from lidar. User must have grid.file_object attribute
        set first using the read_lidar method.
        """
        self.NoDATA = 0
        self.array = None
        self.grid_type = 'intensity'
        
        las = self.file_object

        # extract minimum and maximum from input file
        self.min = [floor(las.header.min[0]), floor(las.header.min[1])]
        # self.max = [ceil(las.header.max[0]), ceil(las.header.max[1])]
        self.max = [floor(las.header.max[0])+1, floor(las.header.max[1])+1]

        # extract x, y values from file to be used later for density calculations
        las_x, las_y, las_i = las.x, las.y, las.intensity

        # get x axis distance (m) from las file
        dist_x = self.max[0] - self.min[0]

        # get y axis distance (m) from las file
        dist_y = self.max[1] - self.min[1]

        # calculate number of columns for raster grid
        self.col = int(dist_x/self.cell_size)
        # self.col += 1  # add one to avoid rounding issues

        # Apply -1 to have negative y resolution for raster
        ycell = -1 * self.cell_size

        # calculate number of rows for raster grid
        self.row = int(dist_y/self.cell_size)
        # self.row += 1  # add one to avoid rounding issues

        # Create empty numpy array to write values to
        count = np.zeros((self.row, self.col)).astype(np.int32)

        # Aggregate intensity values
        int_sum = np.zeros((self.row, self.col)).astype(np.int32)

        # Scale or "project" values  of LiDAR data to grid,
        scale_x = (las_x - self.min[0]) / self.cell_size
        scale_y = (las_y - self.min[1]) / ycell

        # change type to integer and save as variables to use for index values
        ind_x = scale_x.astype(np.int32)
        ind_y = scale_y.astype(np.int32)

        # Loop through LiDAR point records, count, and add to raster grid
        if self.verbose:
            print(f'\nCalculating intensity for {os.path.splitext(os.path.basename(self.path))[0]}...')

        # Runtime bottleneck - This is O(n) WRT the number of points in the point cloud
        # We can't do better than this withouth sacrificing accuracy
        for x, y, i in np.nditer([ind_x, ind_y, las_i]):
            count[y, x] += 1
            int_sum[y, x] += i

        # Fill areas areas lacking data with 1 to avoid divide by zero error
        count_noZero = (np.where(count > 0, count, 1)).astype(np.int32)

        # calculate intensity
        int_avg = (int_sum / count_noZero).astype(np.int32)

        # scale intensity grid to 8bit unsigned integers
        # int_scaled = (int_avg / 256).astype(np.int32)

        # Interpolate 0 values in array to avoid any holes in data
        self.array = np.where(np.logical_and(int_avg > 1, int_avg != 0),
                              int_avg, self.NoDATA).astype(np.int32)

        self.array = np.delete(self.array, 0, axis=0)
        # lazy fix for top row of cells being generated outside AOI inexplicably


    def create_grid(self, last_grid=False, filter_by_class=None, class_grid=False, intensity_grid=False, save_grid=False):
        """
        Method to create grid from lidar. User must have grid.file_object attribute
        set first using the read_lidar method. User can filter by last return, class code,
        or intensity

        filter by keyword argument parameters, either filter by class, last return, both, or neither
        then extract x, y, z values from file to be used later for density calculations
        """

        with self.file_object as las:

# -------------------------------------------------------------------------
# Last Return Density Grid ------------------------------------------------
# -------------------------------------------------------------------------
            if last_grid:

                self.NoDATA_last = 0
                self.last = True

                # extract minimum and maximum from input file; floor min & ceiling max
                self.min = [floor(las.header.min[0]), floor(las.header.min[1])]
                # self.max = [ceil(las.header.max[0]), ceil(las.header.max[1])]
                self.max = [floor(las.header.max[0])+1, floor(las.header.max[1])+1]

                # get x and axis distance (m) from las file
                dist_x = self.max[0] - self.min[0]
                dist_y = self.max[1] - self.min[1]
                # calculate number of columns for raster grid
                self.col = int(dist_x/self.cell_size)
                # self.col += 1
                self.row = int(dist_y/self.cell_size)
                # self.row += 1
                
                # Create empty numpy array to write values to
                count = np.zeros((self.row, self.col)).astype(np.int32)

                # Apply -1 to have negative y resolution for raster
                ycell = -1 * self.cell_size

                las_return_x, las_return_y = (
                    las.x[las.num_returns == las.return_num],
                    las.y[las.num_returns == las.return_num]
                    )
                # Scale or "project" values  of LiDAR data to grid,
                scale_x = (las_return_x - self.min[0]) / self.cell_size
                scale_y = (las_return_y - self.min[1]) / ycell  
        

                # change type to integer and save as variables to use for index values
                ind_x = scale_x.astype(np.int32)
                ind_y = scale_y.astype(np.int32)

                # Loop through LiDAR point records, count, and add to raster grid
                # if self.verbose:
                    # print(f'\nCalculating density for {os.path.splitext(os.path.basename(self.path))[0]}...')

                # Runtime bottleneck - This is O(n) WRT the number of points in the point cloud
                # We can't do better than this withouth sacrificing accuracy
                for x, y in np.nditer([ind_x, ind_y]):
                    count[y, x] += 1

                # Fill areas areas lacking data with keyword argument specified no data value
                count_noData = (np.where(count > 0, count, self.NoDATA_last)).astype(np.int32)

                # calculate density
                self.last_array = (count_noData / self.cell_size).astype(np.int32)
                # lazy fix for top row of cells being generated outside AOI inexplicably
                self.last_array = np.delete(self.last_array, 0, axis=0)


# -------------------------------------------------------------------------
# Ground Density Grid -----------------------------------------------------
# -------------------------------------------------------------------------

            if class_grid:
                
                self.NoDATA_ground = 0
                self.ground = True

                # extract minimum and maximum from input file; floor min & ceiling max
                self.min = [floor(las.header.min[0]), floor(las.header.min[1])]
                # self.max = [ceil(las.header.max[0]), ceil(las.header.max[1])]
                self.max = [floor(las.header.max[0])+1, floor(las.header.max[1])+1]

                # get x and axis distance (m) from las file
                dist_x = self.max[0] - self.min[0]
                dist_y = self.max[1] - self.min[1]
                # calculate number of columns for raster grid
                self.col = int(dist_x/self.cell_size)
                # self.col += 1
                self.row = int(dist_y/self.cell_size)
                # self.row += 1
                
                # Create empty numpy array to write values to
                count = np.zeros((self.row, self.col)).astype(np.int32)

                # Apply -1 to have negative y resolution for raster
                ycell = -1 * self.cell_size
            
                las_class_x, las_class_y = (
                    las.x[las.Classification == filter_by_class],
                    las.y[las.Classification == filter_by_class]
                    )

                scale_x = (las_class_x - self.min[0]) / self.cell_size
                scale_y = (las_class_y - self.min[1]) / ycell

                # change type to integer and save as variables to use for index values
                ind_x = scale_x.astype(np.int32)
                ind_y = scale_y.astype(np.int32)

                # Loop through LiDAR point records, count, and add to raster grid
                # if self.verbose:
                    # print(f'\nCalculating density for {os.path.splitext(os.path.basename(self.path))[0]}...')

                # Runtime bottleneck - This is O(n) WRT the number of points in the point cloud
                # We can't do better than this withouth sacrificing accuracy
                for x, y in np.nditer([ind_x, ind_y]):
                    count[y, x] += 1

                # Fill areas areas lacking data with keyword argument specified no data value
                count_noData = (np.where(count > 0, count, self.NoDATA_ground)).astype(np.int32)

                # calculate density
                self.ground_array = (count_noData / self.cell_size).astype(np.int32)
                # lazy fix for top row of cells being generated outside AOI inexplicably
                self.ground_array = np.delete(self.ground_array, 0, axis=0)

# -------------------------------------------------------------------------
# Intensity Grid ----------------------------------------------------------
# -------------------------------------------------------------------------
            if intensity_grid:

                self.NoDATA_int = 0
                self.intensity = True

                # extract minimum and maximum from input file
                self.min = [floor(las.header.min[0]), floor(las.header.min[1])]
                # self.max = [ceil(las.header.max[0]), ceil(las.header.max[1])]
                self.max = [floor(las.header.max[0])+1, floor(las.header.max[1])+1]

                # extract x, y values from file to be used later for density calculations
                las_x, las_y, las_i = las.x, las.y, las.intensity

                # get x axis distance (m) from las file
                dist_x = self.max[0] - self.min[0]

                # get y axis distance (m) from las file
                dist_y = self.max[1] - self.min[1]

                # calculate number of columns for raster grid
                self.col = int(dist_x/self.cell_size)
                # self.col += 1  # add one to avoid rounding issues

                # Apply -1 to have negative y resolution for raster
                ycell = -1 * self.cell_size

                # calculate number of rows for raster grid
                self.row = int(dist_y/self.cell_size)
                # self.row += 1  # add one to avoid rounding issues

                # Create empty numpy array to write values to
                count = np.zeros((self.row, self.col)).astype(np.int32)

                # Aggregate intensity values
                int_sum = np.zeros((self.row, self.col)).astype(np.int32)

                # Scale or "project" values  of LiDAR data to grid,
                scale_x = (las_x - self.min[0]) / self.cell_size
                scale_y = (las_y - self.min[1]) / ycell

                # change type to integer and save as variables to use for index values
                ind_x = scale_x.astype(np.int32)
                ind_y = scale_y.astype(np.int32)

                # Loop through LiDAR point records, count, and add to raster grid
                # if self.verbose:
                    # print(f'\nCalculating intensity for {os.path.splitext(os.path.basename(self.path))[0]}...')

                # Runtime bottleneck - This is O(n) WRT the number of points in the point cloud
                # We can't do better than this withouth sacrificing accuracy
                for x, y, i in np.nditer([ind_x, ind_y, las_i]):
                    count[y, x] += 1
                    int_sum[y, x] += i

                # Fill areas areas lacking data with 1 to avoid divide by zero error
                count_noZero = (np.where(count > 0, count, 1)).astype(np.int32)

                # calculate intensity
                int_avg = (int_sum / count_noZero).astype(np.int32)

                # scale intensity grid to 8bit unsigned integers
                # int_scaled = (int_avg / 256).astype(np.int32)

                # Interpolate 0 values in array to avoid any holes in data
                self.intensity_array = np.where(np.logical_and(int_avg > 1, int_avg != 0),
                                    int_avg, self.NoDATA_int).astype(np.int32)

                self.intensity_array = np.delete(self.intensity_array, 0, axis=0)
                # lazy fix for top row of cells being generated outside AOI inexplicably


    def save_grid(self, gridtype='unspecified', outdir=None, output_format='tif'):
        """Method to save grid to output in either Geotiff or ESRI ASCII format"""
        tgt = f'{outdir}\\{gridtype}_grid.{output_format}'
        
        if 'intensity' in gridtype:
            grid = self.intensity_array
            self.NoDATA = self.NoDATA_int
        elif 'ground' in gridtype:
            grid = self.ground_array
            self.NoDATA = self.NoDATA_ground
        elif 'last' in gridtype:
            grid = self.last_array
            self.NoDATA = self.NoDATA_last

        if output_format == 'asc':
            # Create ASCII header
            header = "ncols {}\n".format(grid.shape[1])
            header += "nrows {}\n".format(grid.shape[0])
            header += "xllcorner {}\n".format(floor(self.min[0]))
            header += "yllcorner {}\n".format(ceil(self.min[1]))
            header += "cellsize {}\n".format(self.cell_size)
            header += "NODATA_value {}\n".format(self.NoDATA)

            # Open output file, write header, save array to output file
            if grid is not None:
                with open(tgt, "wb") as f:
                    f.write(bytes(header, 'UTF-8'))
                    # fmt string to output integers only
                    np.savetxt(f, grid, fmt="%.0f")
            else:
                raise Exception('No array attribute found for grid object')

        elif output_format == 'tif':
            # create new raster and write array to image
            driver = gdal.GetDriverByName('GTiff')  # get geotiff driver
            if grid is not None:
                out_img = driver.Create( # create raster object
                    tgt, self.col, self.row, 1, gdal.GDT_Int32,
                    options=['COMPRESS=LZW', 'NUM_THREADS=ALL_CPUS'])
                out_img.SetGeoTransform(( # set positional parameters
                    floor(self.min[0]), self.cell_size, 0,
                    ceil(self.max[1]), 0, self.cell_size*-1))
                out_band = out_img.GetRasterBand(1) # get band from raster object
                out_band.SetNoDataValue(self.NoDATA)
                out_band.WriteArray(grid) # write array to raster band
                # write wkt proj metadata to image if it exists
                if self.crs is not None:
                    out_img.SetProjection(self.crs)
                out_band.FlushCache()
            else:
                raise Exception('No array attribute found for grid object')


if __name__ == '__main__':
    main()
import ogr
import os

shapefile = (
    r"D:\test\las_test_files\classified_and_tiled\void_test\fake_AOI.shp"
)  # test file
breaks = (
    r"D:\test\las_test_files\classified_and_tiled\void_test\fake_breaklines.shp"
)  # test file

"""
This program creates a mask made from a cutout of the AOI within
the bounding box of the AOI, and the breaklines. It is intended to be
used with void detection tool within the LiDAR Quality Control Suite
(LiQCs).

Input:      > breaklines for a project in shapefile format
            > The AOI of a project in shapefile format

Output:     > Mask for project in shapefile format

Written by: Jordan Godau
"""


def create_mask(breaklines=None, project_area=None, outdir=None):

    driver = ogr.GetDriverByName("ESRI Shapefile")  # Create driver object
    inds = driver.Open(project_area, 0)  # open shapefile in read mode
    inLayer = inds.GetLayer()  # call get layer method (mandatory)
    inBreaks = driver.Open(breaklines, 0)  # open shapefile in read mode
    breaksLayer = inBreaks.GetLayer()  # call get layer method (mandatory)
    out = os.path.join(outdir, "mask.shp")  # define output name

    # steal prj file from input and write new one for output
    prj_file = f"{os.path.splitext(project_area)[0]}.prj"  # define prj file out name
    with open(prj_file, "r") as f:  # context manager
        contents = f.read()  # read contents
        prj = open(os.path.splitext(out)[0] + ".prj", "w+")  # open new prj file
        prj.write(contents)  # write contents to new file
        prj.close()  # close new file

    # collect bounding geometry
    bounds_list = []  # empty list for bounds
    geom_list = []  # empty list for feature geometry
    for feature in inLayer:  # for geometry feature in layer
        geom = feature.GetGeometryRef()  # not sure what this does yet
        geom_list.append(geom.GetGeometryRef(0).GetPoints())
        bounds_list.append(geom.GetEnvelope())  # get bounds of each geometry feature

    # Create multipolygon of all bounding boxes with cut outs of the geometry they bound
    multipolygon = ogr.Geometry(ogr.wkbMultiPolygon)  # define empty multipolygon
    outRing = ogr.Geometry(ogr.wkbLinearRing)  # Create geometry ring object
    inRing = ogr.Geometry(ogr.wkbLinearRing)  # define inner ring object
    poly = ogr.Geometry(ogr.wkbPolygon)  # create polygon object
    for b, g in zip(bounds_list, geom_list):  # for bounds of each geometry feature
        outRing.AddPoint(b[0] - 2300, b[2] - 2300)  # add points with 2300m buffer
        outRing.AddPoint(b[0] - 2300, b[3] + 2300)  # add points with 2300m buffer
        outRing.AddPoint(b[1] + 2300, b[3] + 2300)  # add points with 2300m buffer
        outRing.AddPoint(b[1] + 2300, b[2] - 2300)  # add points with 2300m buffer
        outRing.AddPoint(b[0] - 2300, b[2] - 2300)  # close polygon on beginning point
        for p in g:  # for vertices in geometry
            inRing.AddPoint(p[0], p[1])  # add vertex to inner ring
        poly.AddGeometry(outRing)  # Polygon inherits outer ring geometry
        poly.AddGeometry(inRing)  # Polygon inherits inner ring geometry
        multipolygon.AddGeometry(poly)  # add polygon to multipolygon

    # get vertices for breaklines
    breakPoints = []
    for feature in breaksLayer:  # for geometry feature in layer
        geomBreaks = feature.GetGeometryRef()  # not sure what this does yet
        breakPoints.append(geomBreaks.GetGeometryRef(0).GetPoints())

    # make polygon feature from breaklines
    multipolyBreaks = ogr.Geometry(ogr.wkbMultiPolygon)  # define empty multipolygon
    ringBreaks = ogr.Geometry(ogr.wkbLinearRing)  # define inner ring object
    polyBreaks = ogr.Geometry(ogr.wkbPolygon)  # create polygon object
    for t in breakPoints:
        for p in t:
            ringBreaks.AddPoint(p[0], p[1])
        polyBreaks.AddGeometry(ringBreaks)
        multipolyBreaks.AddGeometry(polyBreaks)

    # calculate union for breaks and bounds with AOI hole cutout
    union = multipolygon.Union(multipolyBreaks)

    # save bounds
    outDriver = ogr.GetDriverByName("ESRI Shapefile")

    # Remove output shapefile if it already exists
    if os.path.exists(out):
        outDriver.DeleteDataSource(out)

    # create the datasource for the output shapefile
    outDataSource = outDriver.CreateDataSource(out)  # create datasource object
    outLayer = outDataSource.CreateLayer(
        "mask", geom_type=ogr.wkbMultiPolygon
    )  # create new layer

    # add ID field
    idField = ogr.FieldDefn("id", ogr.OFTInteger)  # define field
    outLayer.CreateField(idField)  # create field

    # create feature and set values
    featureDefn = (
        outLayer.GetLayerDefn()
    )  # mandatory, don't know what this actually does
    feature = ogr.Feature(featureDefn)  # same as above
    feature.SetGeometry(union)  # set geometry for feature class
    feature.SetField("id", 1)  # set field for feature class
    outLayer.CreateFeature(feature)  # create new feature
    feature = None  # remove feature from memory

    # Save and close DataSource
    inds = None
    inBreaks = None
    outDataSource = None


def main():
    create_mask(breaklines=breaks, project_area=shapefile, outdir="D:\\")


if __name__ == "__main__":
    main()


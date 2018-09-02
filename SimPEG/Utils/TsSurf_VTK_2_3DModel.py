#%%
from SimPEG import Mesh
import numpy as np
import time as tm
import vtk, vtk.util.numpy_support as npsup
import re

"""
Loads in a triangulated surface from Gocad (*.ts) and use VTK to transfer onto
a 3D mesh.

New scripts to be added to basecode
Require VTK installed:

For Python 3 users:
conda install -c clinicalgraphics vtk

"""


work_dir = 'C:/Users/DominiqueFournier/Dropbox/UBC/ResearchDomFournier/Paper_II_StructuralConstraints/Notebook/'

mshfile = 'Mesh_20m.msh'
outFile = 'TopoFold_20m.den'
# Load mesh file
mesh = Mesh.TensorMesh.readUBC(work_dir + mshfile)

# Load in observation file
#[B,M,dobs] = PF.BaseMag.readUBCmagObs(obsfile)

# Read in topo surface
topsurf = work_dir + "Topo.ts" #work_dir+'CDED_Lake_Coarse.ts'
geosurf = [[work_dir+'Fold_Tilted_noProp_top.ts', True, False],
           [work_dir+'Fold_Tilted_noProp.ts', True, False]]

# Background density
bkgr = 0
airc = -100

# Units
vals = np.asarray([0.15, 0, 3, 4, 5, 6, 7])

#%% Script starts here
# # Create a grid of observations and offset the z from topo

model = np.ones(mesh.nC) * bkgr
# Load GOCAD surf
#[vrtx, trgl] = PF.BaseMag.read_GOCAD_ts(tsfile)
# Find active cells from surface


def read_GOCAD_ts(tsfile):
    """Read GOCAD triangulated surface (*.ts) file
    INPUT:
    tsfile: Triangulated surface

    OUTPUT:
    vrts : Array of vertices in XYZ coordinates [n x 3]
    trgl : Array of index for triangles [m x 3]. The order of the vertices
            is important and describes the normal
            n = cross( (P2 - P1 ) , (P3 - P1) )


    Created on Jan 13th, 2016

    Author: @fourndo
    """

    fid = open(tsfile, 'r')
    line = fid.readline()

    # Skip all the lines until the vertices
    while re.match('TFACE', line) == None:
        line = fid.readline()

    line = fid.readline()
    vrtx = []

    # Run down all the vertices and save in array
    while re.match('VRTX', line):
        l_input = re.split('[\s*]', line)
        temp = np.array(l_input[2:5])
        vrtx.append(temp.astype(np.float))

        # Read next line
        line = fid.readline()

    vrtx = np.asarray(vrtx)

    # Skip lines to the triangles
    while re.match('TRGL', line) == None:
        line = fid.readline()

    # Run down the list of triangles
    trgl = []

    # Run down all the vertices and save in array
    while re.match('TRGL', line):
        l_input = re.split('[\s*]', line)
        temp = np.array(l_input[1:4])
        trgl.append(temp.astype(np.int))

        # Read next line
        line = fid.readline()

    trgl = np.asarray(trgl)

    return vrtx, trgl


def gocad2vtk(gcFile, mesh, bcflag, inflag):
    """"
    Function to read gocad polystructure file and output indexes of
    mesh with in the structure.

    """
    print("Reading GOCAD ts file...")
    vrtx, trgl = read_GOCAD_ts(gcFile)
    # Adjust the index
    trgl = trgl - 1

    # Make vtk pts
    ptsvtk = vtk.vtkPoints()
    ptsvtk.SetData(npsup.numpy_to_vtk(vrtx, deep=1))

    # Make the polygon connection
    polys = vtk.vtkCellArray()
    for face in trgl:
        poly = vtk.vtkPolygon()
        poly.GetPointIds().SetNumberOfIds(len(face))
        for nrv, vert in enumerate(face):
            poly.GetPointIds().SetId(nrv, vert)
        polys.InsertNextCell(poly)

    # Make the polydata, structure of connections and vrtx
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(ptsvtk)
    polyData.SetPolys(polys)

    # Make implicit func
    ImpDistFunc = vtk.vtkImplicitPolyDataDistance()
    ImpDistFunc.SetInput(polyData)

    # Convert the mesh
    vtkMesh = vtk.vtkRectilinearGrid()
    vtkMesh.SetDimensions(mesh.nNx, mesh.nNy, mesh.nNz)
    vtkMesh.SetXCoordinates(npsup.numpy_to_vtk(mesh.vectorNx, deep=1))
    vtkMesh.SetYCoordinates(npsup.numpy_to_vtk(mesh.vectorNy, deep=1))
    vtkMesh.SetZCoordinates(npsup.numpy_to_vtk(mesh.vectorNz, deep=1))
    # Add indexes
    vtkInd = npsup.numpy_to_vtk(np.arange(mesh.nC), deep=1)
    vtkInd.SetName('Index')
    vtkMesh.GetCellData().AddArray(vtkInd)

    extractImpDistRectGridFilt = vtk.vtkExtractGeometry()
    extractImpDistRectGridFilt.SetImplicitFunction(ImpDistFunc)
    extractImpDistRectGridFilt.SetInputData(vtkMesh)

    if bcflag is True:
        extractImpDistRectGridFilt.ExtractBoundaryCellsOn()

    else:
        extractImpDistRectGridFilt.ExtractBoundaryCellsOff()

    if inflag is True:
        extractImpDistRectGridFilt.ExtractInsideOn()

    else:
        extractImpDistRectGridFilt.ExtractInsideOff()

    print("Extracting indices from grid...")
    # Executing the pipe
    extractImpDistRectGridFilt.Update()

    # Get index inside
    insideGrid = extractImpDistRectGridFilt.GetOutput()
    insideGrid = npsup.vtk_to_numpy(insideGrid.GetCellData().GetArray('Index'))

    # Return the indexes inside
    return insideGrid

# Loop through the surfaces and assign value in model
for ii in range(len(geosurf)):
    tin = tm.time()
    print("Computing indices with VTK: " + geosurf[ii][0])
    indx = gocad2vtk(geosurf[ii][0], mesh,
                     bcflag=geosurf[ii][1], inflag=geosurf[ii][2])
    print("VTK operation completed in " + str(tm.time() - tin) + " sec")

    model[indx] = vals[ii]

# Last step to remove aircells
if topsurf:
    indx = gocad2vtk(topsurf, mesh, bcflag=False, inflag=True)
    actv = np.zeros(mesh.nC)
    actv[indx] = 1

    model[actv == 1] = airc

Mesh.TensorMesh.writeModelUBC(mesh, work_dir+outFile, model)

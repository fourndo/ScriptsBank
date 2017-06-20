clear all
dsep = '\';
work_dir = 'C:\Users\DominiqueFournier\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\1D\DIGHEM_TKC_Line10';
data1D = 'DIGHEM_line10.obs';
topofile = 'CDED_076c05_NAD27.topo';
meshfile = '..\UBC_mesh_small.msh';

% Create project and load core files
proj = GIFproject;
proj.setWorkDir(work_dir);
mesh = mesh3D(proj); mesh.readFile([proj.workDir dsep meshfile]);
topo = TOPOdata(proj); topo.readFile([proj.workDir dsep topofile]);

% Create EM1DFM inversion object and pair
inv = EM1DFMinversion(proj);
data = inv.readDobs([proj.workDir dsep data1D]);
inv.setMesh(mesh);
inv.setTopo(topo);




from gridData import Grid
import pandas as pd
from datatable import fread
from scipy.spatial import cKDTree
import networkx as nx
from biopandas.pdb import PandasPdb
import argparse
from pandas.api.types import is_numeric_dtype
##参数
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dx", help="Path of .dx file")
parser.add_argument("-s", "--dms", help="Path of dms output file")
parser.add_argument("-p", "--pdb", help="Path of .pdb file")
parser.add_argument("-o", "--out", help="Path of output file")
args = parser.parse_args()
pdb_file = args.pdb
dms_file = args.dms
dx_file = args.dx
out_file = args.out

##格点文件
g = Grid(dx_file)
data =  g.grid
grid_spacing = 1
###起始点
with open(dx_file, 'r') as f:
    lines = f.readlines()
grid_origin = [float(i) for i in lines[5].split()[-3:]]

x_all = []
y_all = []
z_all = []
all_value = []
for i in range(0,97):
    for j in range(0,97):
        for k in range(0,97):
            value = data[i, j, k]
            x = grid_origin[0] + i*grid_spacing
            y = grid_origin[1] + j*grid_spacing
            z = grid_origin[2] + k*grid_spacing
            x_all.extend([x])
            y_all.extend([y])
            z_all.extend([z])
            all_value.extend([value])

##筛选电势
df = pd.DataFrame(list(zip(x_all, y_all, z_all, all_value)), columns =['x', 'y', 'z', 'potential'])
df = df[df.potential>2].reset_index(drop=True)
##表面文件
dms = fread(dms_file,sep = " ",fill=True).to_pandas()
if not is_numeric_dtype(dms["C5"]):
    dms["C51"] = pd.to_numeric(dms["C5"], errors='coerce')
    dms_c5_nan = dms[dms['C51'].isna()]
    dms_drop_nan = dms.dropna(subset=['C51'])
    dms_c5_nan[['C21', 'C31']] = dms_c5_nan['C2'].str.split('-', expand=True)
    dms_c5_nan['C31'] = pd.to_numeric(dms_c5_nan['C31'])
    dms_c5_nan["C32"] = -1 * dms_c5_nan.C31
    dms_c5_nan = dms_c5_nan[["C0","C1","C21","C32","C3","C4","C5"]]
    dms_c5_nan.columns = ["C0","C1","C2","C3","C4","C5","C6"]
    dms_drop_nan = dms_drop_nan[["C0","C1","C2","C3","C4","C51","C6"]]
    dms_drop_nan.columns = ["C0","C1","C2","C3","C4","C5","C6"]
    dms = pd.concat([dms_drop_nan,dms_c5_nan], ignore_index=True)
dms = dms[dms.C6 != "A"]
###表面的符合电势的点
tree = cKDTree(dms[['C3', 'C4', 'C5']])
indices = tree.query_ball_point(df[['x', 'y', 'z']], r=1)
need_point = [i for i in range(len(indices)) if len(indices[i])>0]
df_need = df.iloc[need_point].reset_index(drop=True)
###找patches
tree1 = cKDTree(df_need[['x', 'y', 'z']])
indices1 = tree1.query_ball_point(df_need[['x', 'y', 'z']], r=1)
##转化为 dict of list
indices1_edge = {}
for i in range(len(df_need)):
    indices1_edge[i] = indices1[i]
    
G = nx.from_dict_of_lists(indices1_edge)
patches = list(nx.connected_components(G))
pathches_size = [len(i) for i in patches]
top_patches_index = [pathches_size.index(i) for i in sorted(pathches_size, reverse=True)[0:3]]
##前3个patch
patch1 = df_need.iloc[list(patches[top_patches_index[0]])].reset_index(drop=True)
patch2 = df_need.iloc[list(patches[top_patches_index[1]])].reset_index(drop=True)
patch3 = df_need.iloc[list(patches[top_patches_index[2]])].reset_index(drop=True)
##pdb
egfr_pdb = PandasPdb().read_pdb(pdb_file)
egfr_pdb = egfr_pdb.df['ATOM']

tree_patch1 = cKDTree(patch1[['x', 'y', 'z']])
indices_patch1 = tree_patch1.query_ball_point(egfr_pdb[['x_coord', 'y_coord', 'z_coord']], r=1)
patch1_aa = egfr_pdb.iloc[[i for i in range(len(indices_patch1)) if len(indices_patch1[i])>0]][["residue_name","chain_id","residue_number"]].drop_duplicates()

tree_patch2 = cKDTree(patch2[['x', 'y', 'z']])
indices_patch2 = tree_patch2.query_ball_point(egfr_pdb[['x_coord', 'y_coord', 'z_coord']], r=1)
patch2_aa = egfr_pdb.iloc[[i for i in range(len(indices_patch2)) if len(indices_patch2[i])>0]][["residue_name","chain_id","residue_number"]].drop_duplicates()

tree_patch3 = cKDTree(patch3[['x', 'y', 'z']])
indices_patch3 = tree_patch3.query_ball_point(egfr_pdb[['x_coord', 'y_coord', 'z_coord']], r=1)
patch3_aa = egfr_pdb.iloc[[i for i in range(len(indices_patch3)) if len(indices_patch3[i])>0]][["residue_name","chain_id","residue_number"]].drop_duplicates()

patch1_aa["patch"] = "p1"
patch2_aa["patch"] = "p2"
patch3_aa["patch"] = "p3"

all_res = pd.concat([patch1_aa, patch2_aa, patch3_aa])
all_res.to_csv(out_file, index=False)

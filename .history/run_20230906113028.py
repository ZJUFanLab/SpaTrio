import spatrio
import pandas as pd
import argparse


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser(description='Spatrio argparse')

# load_data
parser.add_argument('--input_path', default=None, help='Path to read in input data')
parser.add_argument('--ref_path', default=None, help='Path to read in reference data')
# process_input
parser.add_argument('--marker_use', default=True, type=str2bool)
parser.add_argument('--top_marker_num', default=100, type=int)
parser.add_argument('--hvg_use', default=False, type=str2bool)
# ot_alignment
parser.add_argument('--alpha', default=0.1, type=float)
parser.add_argument('--dissimilarity', default='scaled_euc', type=str)
parser.add_argument('--k', default=10, type=int)
parser.add_argument('--graph_mode', default='connectivity', type=str)
parser.add_argument('--aware_spatial', default=True, type=str2bool)
parser.add_argument('--aware_multi', default=True, type=str2bool)
parser.add_argument('--aware_power', default=2, type=int)
# assign_coord
parser.add_argument('--top_num', default=None, type=int)
parser.add_argument('--random', default=False, type=str2bool)
# save output
parser.add_argument('--output_path', default='', help='Path to save output data')

args = parser.parse_args()


path = args.input_path
ref_path = args.ref_path
top_num = args.top_num

print('\n','***Spatrio is running***','\n')

print('\n','***STEP1 LOAD DATA***','\n')
print('Loading data...')
spot_ann = spatrio.load_data(path+'/spatial_rna.csv' )
single_ann = spatrio.load_data( path+'/multi_rna.csv' )
spot_meta =  pd.read_csv(path+'/spatial_meta.csv', index_col=0)
spot_meta['type'] = spot_meta['type'].apply(lambda x: str(x))
single_meta =  pd.read_csv(path+'/multi_meta.csv', index_col=0)
single_meta['type'] = single_meta['type'].apply(lambda x: str(x))
spot_meta[['sample']] = 'spot'
single_meta[['sample']] = 'single'
spot_ann.obs['type'] = spot_meta['type']
spot_ann.obs['type'] = spot_ann.obs['type'].astype(object)
single_ann.obs['type'] = single_meta['type']
single_ann.obs['type'] = single_ann.obs['type'].astype(object)
print('Done!')

pos = pd.read_csv(path+'/pos.csv', index_col=0)
emb = pd.read_csv(path+'/emb.csv', index_col=0)
spot_ann.obsm['spatial'] = pos
single_ann.obsm['reduction'] = emb

if ref_path is not None:
    ref_counts = pd.read_csv(ref_path,index_col=0)
    ref_ratios = ref_counts.div(ref_counts.sum(axis=1), axis=0)
    expected_num = pd.DataFrame({'cell_num':ref_counts.sum(axis=1).tolist()},index = ref_counts.index.values)
elif top_num is not None:
    tmp = pd.DataFrame({'spot': spot_ann.obs.index.tolist(), 'spot_type': spot_ann.obs['type']})
    ref_counts = tmp.groupby(['spot', 'spot_type']).size().unstack(fill_value=0)*top_num
    ref_ratios = ref_counts.div(ref_counts.sum(axis=1), axis=0)
    expected_num = None

print('\n','***STEP2 PROCESS DATA***','\n')
print('Processing data...')
data1,data2 = spatrio.process_input(spot_ann,single_ann,
                                    marker_use = args.marker_use,
                                    top_marker_num = args.top_marker_num,
                                    hvg_use = args.hvg_use)
print('Done!')

print('\n','***STEP3 OT ALIGNMENT***','\n')
if expected_num is not None:
    expected_num = expected_num.loc[data1.obs_names]
    tmp_num = expected_num['cell_num'].values
    p_distribution  = tmp_num/sum(tmp_num)
else:
    p_distribution = None
spatrio_decon = spatrio.ot_alignment(adata1 = data1, 
                                     adata2 = data2, 
                                     dissimilarity = args.dissimilarity,
                                     alpha = args.alpha, 
                                     k = args.k,
                                     graph_mode = args.graph_mode,
                                     aware_spatial = args.aware_spatial,
                                     aware_multi = args.aware_multi,
                                     aware_power = args.aware_power,
                                     p_distribution = p_distribution)
print('\n','***STEP4 ASSIGN COORD***','\n')
spatrio_map = spatrio.assign_coord(adata1 = data1,
                                   adata2 = data2,
                                   out_data = spatrio_decon,
                                   top_num = args.top_num,
                                   random = args.random,
                                   expected_num = expected_num)

print('\n','The output was saved as '+str(args.output_path)+'/output.csv')
spatrio_map.to_csv(args.output_path+'/'+'output.csv')

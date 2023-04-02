import spatrio
import numpy as np
import pandas as pd

path = '/home/yph/spatrio_data/simulation/sim_data/snare-seq/circle'

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

pos = pd.read_csv(path+'/pos.csv', index_col=0)
emb = pd.read_csv(path+'/emb.csv', index_col=0)
spot_ann.obsm['spatial'] = pos
single_ann.obsm['reduction'] = emb

test_num = 10
myalpha = [0,0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
pseudocount=  range(0,6)
accuracy_path = path+'/accuracy_df_aware.csv'
pearson_path = path+'/pearson_df_aware.csv'

alpha_list = []
for i in myalpha:
    alpha_list.append(("alpha_"+str(i)))
accuracy_df_all = pd.DataFrame(columns=alpha_list)
pearson_df_all = pd.DataFrame(columns=alpha_list)

for pc_num in pseudocount:
    for j in range(1,test_num+1):
        accuracy = []
        pearson = []
        data1,data2 = spatrio.process_input(spot_ann,single_ann,marker_use=False,min_cells=0)
        data1 = spatrio.simulate_gene_exp(data1, pc = pc_num, factor= 1)
        for i in myalpha:
            print("### pseudocount = "+str(pc_num)+' rep'+str(j)+' alpha = '+str(i))
            spatrio_decon = spatrio.ot_alignment(adata1 = data1, 
                                                 adata2 = data2, 
                                                 alpha = i, 
                                                 aware_power = 2,
                                                 aware_spatial = True,
                                                 aware_multi=True,
                                                 use_gpu=False)
            spatrio_map = spatrio.assign_coord(adata1 = data1,
                                               adata2 = data2,
                                               out_data = spatrio_decon,
                                               random=False,
                                               top_num=5)

            accuracy_s = spatrio.accuracy_cal(spatrio_map)
            pearson_s = spatrio.pearson_cal(spatrio_map, spot_ann, single_ann)

            accuracy.append(accuracy_s)
            pearson.append(pearson_s)

        accuracy_df = pd.DataFrame(accuracy,index=alpha_list,columns=["rep_"+str(j)]).T
        accuracy_df_all = pd.concat([accuracy_df_all,accuracy_df])
        pearson_df = pd.DataFrame(pearson,index=alpha_list,columns=["rep_"+str(j)]).T
        pearson_df_all = pd.concat([pearson_df_all,pearson_df])

accuracy_df_all['pseudocount'] = np.hstack([[aa]*test_num for aa in pseudocount])
pearson_df_all['pseudocount'] = np.hstack([[aa]*test_num for aa in pseudocount])
accuracy_df_all.to_csv(accuracy_path)
pearson_df_all.to_csv(pearson_path)


###########################
###########################
###########################

test_num = 10
myalpha = [0,0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
pseudocount=  range(0,6)
accuracy_path = path+'/accuracy_df.csv'
pearson_path = path+'/pearson_df.csv'

alpha_list = []
for i in myalpha:
    alpha_list.append(("alpha_"+str(i)))
accuracy_df_all = pd.DataFrame(columns=alpha_list)
pearson_df_all = pd.DataFrame(columns=alpha_list)

for pc_num in pseudocount:
    for j in range(1,test_num+1):
        accuracy = []
        pearson = []
        data1,data2 = spatrio.process_input(spot_ann,single_ann,marker_use=False,min_cells=0)
        data1 = spatrio.simulate_gene_exp(data1, pc = pc_num, factor= 1)
        for i in myalpha:
            print("### pseudocount = "+str(pc_num)+' rep'+str(j)+' alpha = '+str(i))
            spatrio_decon = spatrio.ot_alignment(adata1 = data1, 
                                                 adata2 = data2, 
                                                 alpha = i, 
                                                 aware_power = 2,
                                                 aware_spatial = False,
                                                 aware_multi=False,
                                                 use_gpu=False)
            spatrio_map = spatrio.assign_coord(adata1 = data1,
                                               adata2 = data2,
                                               out_data = spatrio_decon,
                                               random=False,
                                               top_num=5)

            accuracy_s = spatrio.accuracy_cal(spatrio_map)
            pearson_s = spatrio.pearson_cal(spatrio_map, spot_ann, single_ann)

            accuracy.append(accuracy_s)
            pearson.append(pearson_s)

        accuracy_df = pd.DataFrame(accuracy,index=alpha_list,columns=["rep_"+str(j)]).T
        accuracy_df_all = pd.concat([accuracy_df_all,accuracy_df])
        pearson_df = pd.DataFrame(pearson,index=alpha_list,columns=["rep_"+str(j)]).T
        pearson_df_all = pd.concat([pearson_df_all,pearson_df])

accuracy_df_all['pseudocount'] = np.hstack([[aa]*test_num for aa in pseudocount])
pearson_df_all['pseudocount'] = np.hstack([[aa]*test_num for aa in pseudocount])
accuracy_df_all.to_csv(accuracy_path)
pearson_df_all.to_csv(pearson_path)
import streamlit as st
import streamlit.components as stc
import base64 
import time
timestr = time.strftime("%Y%m%d-%H%M%S")
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Make a table to represent a 384 well plate
ls_384well = []

# Number of rows and columns in a 384 well plate
num_plate_rows = 16
num_plate_cols = 24

# Potential row identifiers
alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', \
            'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', \
            'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', \
            'Y', 'Z']
alphabet = alphabet[:num_plate_rows]
# Column identifiers
col_ids = list(range(1, 25))

well_dict = {}
plate_384well_ls = []

# Fill in the well_dict
for row in range(num_plate_rows):
    # Get current row identifier from the alphabet
    curr_row_id = alphabet[row]
    row_ls = []
    
    # Get all ids in the row
    for curr_col_id in col_ids:
        if curr_col_id < 10:
            row_ls.append(curr_row_id + '0' + str(curr_col_id))
        else:
            row_ls.append(curr_row_id + str(curr_col_id))
    
    # Retain well ids for the row in the current row identifier
    well_dict[curr_row_id] = row_ls
    plate_384well_ls += row_ls
    
plate_384_format = pd.DataFrame.from_dict(well_dict, orient = 'index', columns = col_ids)

# From the pooled hits, identify wells from source plate
def pool_to_source(hits, pooling_type = 'plate'):
    # Determine how plates were pooled
    p_type = pooling_type.lower()
    
    hit_ls = []
    
    # Plates were pooled by quadrant -- Plate 1's A01, A02, B01, and B02 were pooled
    if 'q' in p_type:
        for ind, hit in hits.iterrows():
            pooled_plate_num = hit['SourcePlate']
            pooled_well = hit['SourceWell']
            pooled_row = pooled_well[0]
            pooled_col = int(pooled_well[1:])
            
            # if in row A, C, E, G, I, K, M, O
            #    either plate 1 or 2
            if pooled_row in alphabet[::2]:
                # Get rows in source plate to draw from
                source_rows = [pooled_row, alphabet[alphabet.index(pooled_row) + 1]]
                # if in odd col, then in plate 1
                # if in even col, then in plate 2
                if pooled_col % 2 == 0:
                    source_plate_num = (pooled_plate_num * 4) - (4-2)
                    source_cols = [pooled_col - 1, pooled_col]
                else:
                    source_plate_num = (pooled_plate_num * 4) - (4-1)
                    source_cols = [pooled_col, pooled_col + 1]
            
            # if in row B, D, F, H, J, L, N, P
            # either plate 3 or 4
            else:
                # Get rows in source plate to draw from
                source_rows = [alphabet[alphabet.index(pooled_row) - 1], pooled_row]
                # if in odd col, then in plate 3
                # if in even col, then in plate 4
                if pooled_col % 2 == 0:
                    source_plate_num = (pooled_plate_num * 4) - (4-4)
                    source_cols = [pooled_col - 1, pooled_col]
                else:
                    source_plate_num = (pooled_plate_num * 4) - (4-3)
                    source_cols = [pooled_col, pooled_col + 1]
            
            try:
                for r in source_rows:
                    for c in source_cols:
                        if c < 10:
                            source_well = r + '0' + str(c)
                        else:
                            source_well = r + str(c)
                        hit_ls.append([source_plate_num, source_well])
            except:
                print('Trouble finding source wells from pooled well: %s' % pooled_well)
        
        hit_df = pd.DataFrame(hit_ls, columns = ['SourcePlate', 'SourceWell'])
        return hit_df
                    
    # Plates were pooled by plates -- Plate 1's A01, Plate 2's A01, Plate 3's A01, and Plate 4's A01 were pooled)
    elif 'p' in p_type:
        for ind, hit in hits.iterrows():
            pooled_plate_num = hit['SourcePlate']
            pooled_well = hit['SourceWell']
            
            for i in range(1, 5):
                hit_ls.append([(pooled_plate_num*4) - (4-i), pooled_well])
        hit_df = pd.DataFrame(hit_ls, columns = ['SourcePlate', 'SourceWell'])
        return hit_df
        
    else:
        print('Unknown pool type.')
        return


# Input for deconvolution screen -- hit rearrangement from pooled screen
def deconvolution(hits, pooling_type = 'plate'):
    hits_df = pool_to_source(hits, pooling_type)
    num_rows, num_cols = hits_df.shape
    num_deconv_plates = int(np.ceil(num_rows/384))
    
    # Make deconvolution plate
    deconv_well_ls = []
    deconv_plate_ls = []
    for i in range(num_deconv_plates):
        deconv_well_ls += plate_384well_ls
        curr_plate_ls = [i+1]*384
        deconv_plate_ls += curr_plate_ls
    
    # Make deconvolution plate and well lists the same size as the hits list
    deconv_well_ls = deconv_well_ls[:num_rows]
    deconv_plate_ls = deconv_plate_ls[:num_rows]
    vol_ls = [40]*num_rows
    
    hits_df = hits_df.assign(DestPlate = deconv_plate_ls)
    hits_df = hits_df.assign(DestWell = deconv_well_ls)
    hits_df = hits_df.assign(Vol = vol_ls)
    
    return hits_df

# Input for hit rearrangment -- hit rearrangement from deconvolution screen
def hit_rearrangement(hits):
    hits_df = hits.copy()
    num_rows, num_cols = hits_df.shape
    num_hit_plates = int(np.ceil(num_rows/384))
    
    # Make hit plate
    hit_well_ls = []
    hit_plate_ls = []
    for i in range(num_hit_plates):
        hit_well_ls += plate_384well_ls
        curr_plate_ls = [i+1]*384
        hit_plate_ls += curr_plate_ls
    
    # Make deconvolution plate and well lists the same size as the hits list
    hit_well_ls = hit_well_ls[:num_rows]
    hit_plate_ls = hit_plate_ls[:num_rows]
    vol_ls = [50]*num_rows
    
    hits_df = hits_df.assign(DestPlate = hit_plate_ls)
    hits_df = hits_df.assign(DestWell = hit_well_ls)
    hits_df = hits_df.assign(Vol = vol_ls)
    
    return hits_df

class FileDownloader(object):
	
	def __init__(self, data,filename='myfile',file_ext='txt'):
		super(FileDownloader, self).__init__()
		self.data = data
		self.filename = filename
		self.file_ext = file_ext

	def download(self):
		b64 = base64.b64encode(self.data.encode()).decode()
		new_filename = "{}_{}_.{}".format(self.filename,timestr,self.file_ext)
		href = f'<a href="data:file/{self.file_ext};base64,{b64}" download="{new_filename}">Click Here!!</a>'
		st.markdown(href,unsafe_allow_html=True)

# Following: https://pythonwife.com/file-upload-download-with-streamlit/ to upload and download files
def main():

    st.title('Hamilton programs')

    menu = ['Hit picking', 'Deconvolution', 'Hit rearrangement']

    choice = st.sidebar.selectbox('Select a function:', menu)

    if choice == 'Hit picking':
        st.write('You have chosen: Hit picking')
        st.write('Please load a binding data:')
        
        # If user would like to see an example of how the csv file should look
        plate_ex = [1, 1, 1, 2, 2, 3, 3, 3]
        well_ex = ['A02', ' A06', 'A11', 'B04', 'C06', 'C13', 'A07', 'A21']
        ag_pos_ex = [100, 10000, 50, 1039, 1023, 123, 78, 67]
        ag_neg_ex = [87, 242, 52, 102, 194, 62, 60, 72]
        df_ex = pd.DataFrame([plate_ex, well_ex, ag_pos_ex, ag_neg_ex]).transpose()
        df_ex.columns = ['SourcePlate', 'SourceWell', 'Ag+', 'Ag-']

        if st.checkbox('Show example:'):
            st.write(df_ex)
            
        data_file = st.file_uploader('Upload csv', type = ['csv', 'xlsx'])
        
        if data_file is not None:
            # To see details
            file_details = {'filename': data_file.name, 
                'filetype': data_file.type, 
                'filesize': data_file.size}
        
            data_file_name = data_file.name.split('.')[0]
            if data_file.name.endswith('csv'):
                df = pd.read_csv(data_file)
            elif data_file.name.endswith('xlsx'):
                df = pd.ExcelFile(data_file).parse()
            
            # Plot hit picking
            f, ax = plt.subplots(figsize = (5, 3))
            sns.scatterplot(x = 'Ag+', y = 'Ag-', data = df, ax = ax, s = 100)
            
            df['Fold'] = df['Ag+'].divide(df['Ag-'])
            fold_slider = st.slider('Fold change', 1, 100)
            final_hits = df[df['Fold'] > fold_slider]
            sns.scatterplot(x = 'Ag+', y = 'Ag-', data = final_hits, ax = ax, s = 100, fc = 'red')
            st.pyplot(f)
            st.dataframe(final_hits)
            
            download = FileDownloader(
                final_hits.to_csv(), 
                filename = '%s_hit_list'%data_file_name, 
                file_ext='csv').download()
            
    
    if choice == 'Hit rearrangement':
        st.write('You have chosen: Hit rearrangement')
        st.write('Please load a hit list for rearrangement:')
        
        # If user would like to see an example of how the csv file should look
        plate_ex = [1, 1, 1, 2, 2, 3, 3, 3]
        well_ex = ['A02', ' A06', 'A11', 'B04', 'C06', 'C13', 'A07', 'A21']
        df_ex = pd.DataFrame([plate_ex, well_ex]).transpose()
        df_ex.columns = ['SourcePlate', 'SourceWell']

        if st.checkbox('Show example:'):
            st.write(df_ex)
            
        data_file = st.file_uploader('Upload csv', type = ['csv', 'xlsx'])
        
        if data_file is not None:
            # To see details
            file_details = {'filename': data_file.name, 
                'filetype': data_file.type, 
                'filesize': data_file.size}
        
            data_file_name = data_file.name.split('.')[0]
            if data_file.name.endswith('csv'):
                df = pd.read_csv(data_file)
            elif data_file.name.endswith('xlsx'):
                df = pd.ExcelFile(data_file).parse()
        
            for col in df.columns.values:
                if 'plate' in col.lower():
                    df.rename(columns = {col: 'SourcePlate'},inplace = True)
                elif 'well' in col.lower():
                    df.rename(columns = {col: 'SourceWell'},inplace = True)
                
            hit_df = hit_rearrangement(df)
            st.dataframe(hit_df)
            
            # If resulting file has to be grouped by a specific number of plates for download
            num_plates = st.number_input('Group plates by: ', 0, 12, 0, 1)

            if num_plates == 0:
                st.markdown("#### Download File ###")
                download = FileDownloader(
                    hit_df.to_csv(index = False), 
                    filename = '%s_deconvoluted'%data_file_name, 
                    file_ext='csv').download()

            else:
                total_plates = hit_df['SourcePlate'].max()
                counter = 0
                while counter*num_plates < total_plates:
                    counter += 1
                    df_chunk = hit_df[hit_df['SourcePlate'] <= counter*num_plates]
                    df_chunk = df_chunk[df_chunk['SourcePlate'] > (counter-1)*num_plates]
                    
                    st.markdown("#### Download File for Plates %d-%d ###"%((counter-1)*num_plates+1, counter*num_plates))
                    download = FileDownloader(
                        df_chunk.to_csv(index = False),
                        filename = '%s_deconvoluted_%d'%(data_file_name, counter),
                        file_ext='csv').download()
        
    
    if choice == 'Deconvolution':
        st.write('You have chosen: Deconvolution')
        st.write('Please load a hit list for deconvolution rearrangement:')
        
        # If user would like to see an example of how the csv file should look
        plate_ex = [1, 1, 1, 2, 2, 3, 3, 3]
        well_ex = ['A02', ' A06', 'A11', 'B04', 'C06', 'C13', 'A07', 'A21']
        df_ex = pd.DataFrame([plate_ex, well_ex]).transpose()
        df_ex.columns = ['SourcePlate', 'SourceWell']

        if st.checkbox('Show example:'):
            st.write(df_ex)
            
        data_file = st.file_uploader('Upload csv', type = ['csv', 'xlsx'])
    
        # Select how hybridoma plates were pooled
        pool_choice = st.radio('Plates were pooled by:', ('quadrant', 'plate'))
        
        if data_file is not None:
            # To see details
            file_details = {'filename': data_file.name, 
                'filetype': data_file.type, 
                'filesize': data_file.size}
        
            data_file_name = data_file.name.split('.')[0]
            if data_file.name.endswith('csv'):
                df = pd.read_csv(data_file)
            elif data_file.name.endswith('xlsx'):
                df = pd.ExcelFile(data_file).parse()
                
            for col in df.columns.values:
                if 'plate' in col.lower():
                    df.rename(columns = {col: 'SourcePlate'},inplace = True)
                elif 'well' in col.lower():
                    df.rename(columns = {col: 'SourceWell'},inplace = True)
        
            if pool_choice == 'quadrant':
                deconvoluted_df = deconvolution(df, 'q')
                deconvoluted_df.sort_values(['SourcePlate', 'SourceWell'], inplace = True)
            elif pool_choice == 'plate':
                deconvoluted_df = deconvolution(df, 'p')
            st.dataframe(deconvoluted_df)

            
            # If resulting file has to be grouped by a specific number of plates for download
            num_plates = st.number_input('Group plates by: ', 0, 12, 0, 1)

            if num_plates == 0:
                st.markdown("#### Download File ###")
                download = FileDownloader(
                    deconvoluted_df.to_csv(index = False), 
                    filename = '%s_deconvoluted'%data_file_name, 
                    file_ext='csv').download()

            else:
                total_plates = deconvoluted_df['SourcePlate'].max()
                counter = 0
                while counter*num_plates < total_plates:
                    counter += 1
                    df_chunk = deconvoluted_df[deconvoluted_df['SourcePlate'] <= counter*num_plates]
                    df_chunk = df_chunk[df_chunk['SourcePlate'] > (counter-1)*num_plates]
                    
                    st.markdown("#### Download File for Plates %d-%d ###"%((counter-1)*num_plates+1, counter*num_plates))
                    download = FileDownloader(
                        df_chunk.to_csv(index = False),
                        filename = '%s_deconvoluted_%d'%(data_file_name, counter),
                        file_ext='csv').download()
                    
                    
                

main()

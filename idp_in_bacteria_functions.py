# functions for the IDP in bacteria project

### dependencies ###

from IDP_analysis.packages_import import *


### Read fasta ###

def read_fasta(fasta_file: str) -> pd.DataFrame:
  """Processes raw .fasta files

  Opens a .fasta file, parses the sequences and their IDs into a
  dataframe and returns the dataframe.

  Parameters
  ----------
  fasta_file : str
    the raw .fasta file directory.

  Returns
  -------
  pd.DataFrame
    a dataframe with ID, Sequence and Length columns.

  """

  # open the file
  handle = open(fasta_file, 'r')
  seq_list = list(SeqIO.parse(handle, 'fasta'))
  handle.close()

  # parse data into lists
  ids = [seq_record.id.split('|')[1] for seq_record in seq_list]
  seqs = [str(seq_record.seq) for seq_record in seq_list]
  lens = [len(seq) for seq in seqs]

  # save data into a dataframe
  df = pd.DataFrame({'ID':ids, 'Sequence':seqs, 'Length':lens})

  return df


### Find longest binary IDR ###

def find_longest_binary_IDR(binary_disorder_list: list) -> int:
  """Finds the longest stretch of binary IDR for a
  given protein

  Parameters
  ----------
  binary_disorder_list : list
    a list of binary values indicating whether the
    residue is disordered (1) or not (0)

  Returns
  -------
  int
    a number of disordered residues in the longest
    IDR

  """

  # set max IDR length as zero
  max_IDR_len = 0

  # set current IDR length as zero
  current_IDR_len = 0
  # iterate though each residue
  for disorder in binary_disorder_list:
    if disorder == 1:
      # increase current IDR length if disordered
      current_IDR_len += 1
    else:
      # not part of IDR, save current length if max or continue
      if current_IDR_len > max_IDR_len:
        max_IDR_len = current_IDR_len
      current_IDR_len = 0
  # check last time
  if current_IDR_len > max_IDR_len:
    max_IDR_len = current_IDR_len

  return max_IDR_len


### Align sequences ###

def align_seqs(seqs: list, ids: list, gap_ins_pen: float=-1.0,
               gap_ext_pen: float=-0.5):
  """Multiple sequence alignment using ClustalW

  Parameters
  ----------
  seqs : list
    A list of sequences to be aligned (as str)
  ids : list
    A list of sequence IDs (as str)
  gap_ins_pen : float
    Gap insertion penalty for the alignment.
    -1.0 by default
  gap_ext_pen : float
    Gap extension penalty for the alignment.
    -0.5 by default

  Returns
  -------
  Alignment object

  """

  # create seq records
  records = [SeqRecord(Seq(seq), id) for seq, id in zip(seqs, ids)]
  # save records as fasta
  with open('unaligned.fas', 'w') as handle:
    SeqIO.write(records, handle, 'fasta')
  # align using ClustalW
  cline = ClustalwCommandline(infile='unaligned.fas', outfile='aligned.fas',
                            gapext=gap_ext_pen, gapopen=gap_ins_pen,
                            type='PROTEIN', outorder='INPUT',
                            output='FASTA')
  cline()
  # read the output file
  align = AlignIO.read('aligned.fas', 'fasta')
  return align


### Align disorder ###

def align_disorder(disorder_values: list, seqs: list, ids: list,
                  gap_ins_pen: float=-1.0, gap_ext_pen: float=-0.5) -> tuple:
  """Aligns disorder values based on AA sequences

  Parameters:
  -----------
  disorder_values : list
    A list of disorder values to be plotted
  seqs : list
    A list of sequences to be aligned
  ids : list
    A list of sequence IDs (as str)
  gap_ins_pen : float
    Gap insertion penalty for the alignment.
    -1.0 by default
  gap_ext_pen : float
    Gap extension penalty for the alignment.
    -0.5 by default

  Returns:
  --------
  (sequence_alignment, disorder_alignment)

  """

  aligned_disorders = []
  # align sequences
  alignment = align_seqs(seqs, ids, gap_ins_pen=gap_ins_pen, gap_ext_pen=gap_ext_pen)
  # insert gaps in disorder values
  for a, d in zip(alignment, disorder_values):
    seq_disorder_full = [*d.copy()]
    for i, residue in enumerate(a.seq):
      if residue == '-':
        seq_disorder_full.insert(i, np.nan)
    seq_disorder_full = pd.Series(seq_disorder_full).interpolate().tolist()
    aligned_disorders.append(seq_disorder_full)
  return (alignment, aligned_disorders)


### Get charge / hydrophobicity ###

def get_charge_hydrophobicity(seq: str) -> tuple:
  """Calculates mean absolute charge and hydrophobicity

  Parameters
  ----------
  seq : str
    AA sequence to calculate values for

  Returns:
  --------
  tuple
    tuple with two float values:
    1) mean absolute charge
    2) mean hydrophobicity (scaled 0-1)

  """

  # calculate mean absolute charge
  prot_analysis = ProteinAnalysis(seq)
  mean_abs_charge = abs(prot_analysis.charge_at_pH(7)) / prot_analysis.length

  # get Kyle-Doolitle scale
  KyteDoolitle = gravy_scales.get('KyteDoolitle', -1)
  KyteDoolitle['X'] = 0
  KyteDoolitle['U'] = 0
  # normalize scale to 0-1 range
  scaler = MinMaxScaler()
  KD_norm = {k:v[0] for k, v in zip(KyteDoolitle.keys(), scaler.fit_transform(np.array([*KyteDoolitle.values()]).reshape(-1, 1)))}
  # calculate mean hydrophobicity
  hydrophobicity_sum = sum(KD_norm[aa] for aa in seq)
  mean_hydrophobicity = hydrophobicity_sum / prot_analysis.length

  return mean_abs_charge, mean_hydrophobicity


### Plot aligned cluster ###

def plot_aligned_cluster(df:pd.DataFrame, cluster:str, mav:int=50, errorbar=('ci', 95), scatter:bool=False) -> None:
  """Create a plot of disorder scores for a given cluster using aligned disorder values.
  
  The function generates a plot showing the predicted disorder scores along the sequence positions for the given cluster.
  Each line in the plot represents a species within the cluster, and the plot also includes error bars to indicate the
  standard deviation of disorder scores. The horizontal dashed line at 0.5 serves as a reference point to distinguish
  ordered regions (below 0.5) from disordered regions (above 0.5).

    ----------
    Parameters:
    df (pd.DataFrame) : Dataframe that contains aligned disorder data
    cluster (str) : The identifier of the cluster for which the disorder scores are to be plotted.
    mav (int, optional) : The moving average window size used to smooth the disorder scores. The default value is 50.
    errorbar (tuple or str, optional) : Type of error bar to be used, inherited from sns.lineplot(). Default is ('ci', 95).
    scatter (bool, optional) : Whether to include a scatterplot for each species. Default is False.

    Returns:
    None
    
  """

  # extract disorder data on the specified cluster
  data = df[df['cluster']==cluster]
  aligned_disorders = data['disorder_aligned']
  # calculate the adjusted mav based on data length
  mav = int(np.round(mav * len(aligned_disorders.iloc[0]) / 500)) + 1

  # create figure and axis
  fig, ax = plt.subplots(figsize=(15, 8))
  # set limits and labels
  ax.set_ylim(0, 1)
  ax.set_xlabel('Position in sequence', fontsize=18)
  ax.set_ylabel('Predicted disorder', fontsize=18)
  
  # apply mav
  ys = [np.convolve(ad, np.ones(mav)/mav) for ad in aligned_disorders]
  ys = pd.Series([y[int(mav/2):len(y)-int(mav/2)+1] for y in ys], index=data.index)
  data['ys'] = ys

  # iterate through each group and plot a line with error
  for group in set(data.group):
    group_data = data[data['group']==group]
    color = group_data['color'].iloc[0]
    group_data = pd.DataFrame({group_data['species_tag'].iloc[i]:group_data['ys'].iloc[i] for i in range(len(group_data))})
    group_data['position'] = group_data.index + 1
    group_data = group_data.melt(id_vars=['position'], value_name='disorder', var_name='species_tag')
    sns.lineplot(data=group_data, x='position', y='disorder', errorbar=errorbar, c=color, label=group)
    # scatter for individual species if needed
    if scatter:
      plt.scatter(x=group_data.position, y=group_data.disorder, marker='.', color=color, alpha=0.1)

  # add 0.5 threshold line and legend
  ax.axhline(0.5, ls='--', c='black')
  ax.legend()

  fig.show()


###################################
### OLD FUNCTIONS - no longer used

# ### Read pickle file ###

# def read_pickle_file(file_name: str):
#   """A simple function to read pickle files

#   Parameters
#   ----------
#   file_name : srt
#     location of the pickle file

#   Returns
#   -------
#   any format, but best with pd.DataFrames or lists

#   """

#   print(file_name + ' loading...')
#   with open(file_name, 'rb') as f:
#     result = pickle.load(f)
#   print(file_name + ' loaded!')
#   return result


# ### Binary disorder ###

# def binary_disorder(disorder_list: list, threshold: float=0.5) -> list:
#   """Returns a list of binary values for disordered
#   residues given a threshold

#   Parameters
#   ----------
#   disorder_list : list
#     a list of disorder values
#   threshold : float
#     a binary threshold, 0.5 by default. Disorder score
#     above the threshold means that the residue is
#     considered to be disordered

#   Returns
#   -------
#   list
#     a list of binary values indicating whether the
#     residue is disordered (1) or not (0)

#   """

#   return [1 if dis > threshold else 0 for dis in disorder_list]


# ### Fix -infinitives ###

# def fix_neg_inf(df: pd.DataFrame, replacement=np.nan) -> pd.DataFrame:
#   """Replaces -inf values in a given dataframe

#   Parameters
#   ----------
#   df : pd.DataFrame
#     The dataframe with -inf to be replaced
#   replacement
#     The value to replace -inf with. Default
#     is np.nan

#   Returns:
#   --------
#   pd.DataFrame

#   """

#   df_new = df.copy()
#   df_new[df.astype(float) < 0] = replacement
#   return df_new


# ### Drop unavailable rows ###

# def drop_unavailable(proteome: pd.DataFrame, reset_index: bool=True) -> pd.DataFrame:
#   """Drop rows if not all models are available

#   """

#   mask = [proteome.disorder.iloc[n].shape[1]==4 for n in range(len(proteome))]
#   proteome = proteome[mask]
#   if reset_index == True:
#     proteome.reset_index(inplace=True, drop=True)
#   return proteome


# ### Get additional columns with IDP data ###

# def get_additional_columns(proteome:pd.DataFrame) -> None:
#   """Add Additional Columns to Proteome DataFrame
#   This function calculates additional columns for the proteome DataFrame including:

#   disorder_combined: mean of all disorder models
#   disorder_binary: binary disorder values based on a threshold of 0.5
#   disorder_mean: mean disorder score for each protein
#   disorder_binary_mean: mean binary disorder value (fraction of disorder) for each protein
#   longest_IDR: length of the longest binary intrinsically disordered region for each protein
#   -----------
#   Parameters:
#   proteome (pd.DataFrame): The dataframe containing protein information.

#   Returns:
#   None
#   """

#   # get means of all models
#   proteome['disorder_combined'] = [np.mean(proteome.disorder.iloc[i], axis=1) for i in range(len(proteome))]
#   # get binary disorder values
#   proteome['disorder_binary'] = proteome.disorder_combined.apply(binary_disorder, threshold=0.5)
#   # get mean disorder values for each protein (mean disorder score)
#   proteome['disorder_mean'] = proteome.disorder_combined.apply(np.mean)
#   # get mean binary disorder values for each protein (fraction of disorder)
#   proteome['disorder_binary_mean'] = proteome.disorder_binary.apply(np.mean)
#   # find longest binary IDR for each protein
#   proteome['longest_IDR'] = proteome.disorder_binary.apply(find_longest_binary_IDR)
#   return


# ### Get and process pairwise OMA file ###

# def get_OMA(species_tag:str, common_species_tag:str='THET8', oma_type:str='any') -> pd.DataFrame:
#   """Obtains pairwise OMA and returns a dataframe

#   -----------
#   Parameters:
#   species_tag (str): Tag for species to use in OMA.
#   common_species_tag (str): Tag for the second species to use in OMA. This should be the same
#   for all OMA files. 'THET8' by default.
#   oma_type (str): Type of OMA relationships to keep. 'any' by default. Other possibilities are
#   '1:1', 'n:1', 'm:1' and 'n:m'

#   Returns:
#   pd.DataFrame
#   """

#   # get file
#   path = f'https://omabrowser.org/cgi-bin/gateway.pl?f=PairwiseOrthologs&p1={common_species_tag}&p2={species_tag}&p3=UniProt'
#   try_retrieve_url(path, '_')
#   OMA = pd.read_csv('_', sep='\t')
#   if len(OMA) == 0:
#     path = f'https://omabrowser.org/cgi-bin/gateway.pl?f=PairwiseOrthologs&p1={species_tag}&p2={common_species_tag}&p3=UniProt'
#     try_retrieve_url(path, '_')
#     OMA = pd.read_csv('_', sep='\t')
#   # process file
#   OMA.columns = [common_species_tag, species_tag, 'oma_type', 'OMA_group']
#   if oma_type != 'any':
#     OMA = OMA[OMA['oma_type'] == oma_type]
#   OMA.drop(['OMA_group', 'oma_type'], axis='columns', inplace=True)
#   OMA[common_species_tag] = OMA[common_species_tag].str.replace(common_species_tag, '')
#   OMA[common_species_tag] = OMA[common_species_tag].str.replace('_', '')
#   # OMA = OMA[OMA[common_species_tag].str.len()==6]
#   OMA[species_tag] = OMA[species_tag].str.replace(species_tag, '')
#   OMA[species_tag] = OMA[species_tag].str.replace('_', '')
#   OMA.reset_index(inplace=True, drop=True)
#   return OMA


# ### Read fIDPnn disorder files ###

# def read_fIDPnn_disorder(path, insert_into_df=False, original_df=None):
#   """Processes sequence data from a file and returns DataFrames or updates an existing DataFrame

#     ----------
#     Parameters:
#     path (str) : The path to the input file containing sequence data.
#     insert_into_df (bool, optional) : Whether to insert processed data into an existing DataFrame. Default is False.
#     original_df (pd.DataFrame, optionsl) : Original DataFrame to insert into if insert_into_df=True. Default is None.

#     Returns:
#     None if insert_into_df is True, otherwise a dictionary of DataFrames.
    
#   """
  
#   dfs = {}
#   lists = {}
#   dflines = []

#   with open(path, 'rb') as f:
#     for i, line in enumerate(tqdm(f.readlines())):
#       if line.decode('utf-8')[:8] == 'Sequence':
#         continue
#       elif line.decode('utf-8')[0] == '>':
#         if dflines != []:
#           columns = dflines[0].lower().replace(' ', '_').strip().split(',')
#           df = pd.DataFrame([l.strip().split(',') for l in dflines[1:]], columns=columns)
#           df.iloc[:, 2:] = df.iloc[:, 2:].astype('float')
#           df[df.columns[0]] = df.iloc[:, 0].astype('float')
#           dfs[name] = df
#           lists[name] = df.predicted_score_for_disorder.tolist()
#         name = line.decode('utf-8').split('|')[1]
#         dflines = []
#       else:
#         dflines.append(line.decode('utf-8'))
#   if insert_into_df:
#     m = original_df.disorder_fIDPnn.isna()
#     original_df.loc[m, 'disorder_fIDPnn'] = original_df.ID.map(lists)
#     return
#   else:
#     return dfs



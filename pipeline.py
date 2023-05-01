

import Bio
import Bio.SeqIO
import importlib
import glob
import subprocess

def read_fasta(fname):
    '''
    Read sequences from a fasta file
    fname = file to open
    returns a list with [[filename, name, sequence]]
    '''
    with open(fname, 'r') as f:
        rec = list(Bio.SeqIO.parse(f, 'fasta'))
    df = [[fname, i.name, str(i.seq)] for i in rec]
    col = ['filename','name','sequence']
    df = pd.DataFrame(df, columns=col)
    return df

def trim_adapters():
    # Collect files
    files = glob.glob('data/reads/*.fq*')
    files = [i.split('_')[0] for i in files]
    r1 = [i+'_1_filtered.fq.gz' for i in files]
    r2 = [i+'_2_filtered.fq.gz' for i in files]

    # make an output directory for merged reads
    cmd = 'mkdir trimmed'
    subprocess.run(cmd.split(' '))
    
    # merge the reads
    for i in range(len(files)):
        out1 = 'trimmed/'+r1[i].split('/')[-1]
        out2 = 'trimmed/'+r2[i].split('/')[-1]
        print('trimming',out1)

        # primers found in sample_sheet.csv
        fwd = 'GTGCCAGCMGCCGCGGTAA'
        rev = 'CCGTCAATTCMTTTRAGTTT'
        #rev = 'AAACTYAAAKGAATTGACGG'
        
        #cmd = 'cutadapt --quiet -g '+fwd+' -o '+fout+' '+f
        #subprocess.run(cmd.split(' '))
        #cmd = 'cutadapt --quiet -a '+rev+' -o '+fout+' '+f
        
        cmd = 'cutadapt --quiet -m 100 -g '+fwd+' -G '+rev+' -o '+out1+' -p '+out2+' '+r1[i]+' '+r2[i]
        subprocess.run(cmd.split(' '))

def merge_reads():
    # Collect files
    files = glob.glob('trimmed/*.fq*')
    files = [i.split('_')[0] for i in files]

    # make an output directory for merged reads
    cmd = 'mkdir merged'
    subprocess.run(cmd.split(' '))
    
    # merge the reads
    for f in files:
        # define forward and reverse pairs
        r1 = f+'_1_filtered.fq.gz'
        r2 = f+'_2_filtered.fq.gz'
        fname = f.split('/')[-1]
        fout = 'merged/'+fname+'.fq.gz'
                
        # merge each pair of reads
        print('merging',f)
        cmd = './bin/NGmerge -1 '+r1+' -2 '+r2+' -o '+fout
        subprocess.run(cmd.split(' '))

def cluster_reads():
    # Cluster reads with linclust
    cmd = ['mkdir clusters']
    cmd+= ['mmseqs easy-linclust merged/*.fq.gz clusters/clst tmp']
    #cmd+= ['mmseqs easy-cluster merged/*.fq.gz clusters/clst tmp -c 0.95']
    # we only need clst_rep_seq.fasta
    # remove other files
    cmd+= ['mv clusters/clst_rep_seq.fasta data/clusters.fa']
    cmd+= ['rm clusters/*']
    for c in cmd:
        subprocess.call(c, shell=True)

def parse_PAF(fname):
    '''
    This function parses the information from a PAF file and returns a dataframe
    fname = filename of paf file output from minimap2
    returns a dataframe
    '''
    with open(fname,'r') as f:
        text = f.read().split('\n')

    if len(text[-1])==0:
        text = text[:-1] # removes the last line that is empty    
    data = []

    # information we are parsing
    col1 = ['query_id','q_len','q_start','q_end','orientation',
            'database_id','t_len','t_start','t_end','match','tot','mapq']
    col2 = ['tp','cm','s1','s2','NM','AS','ms','nn','rl','cg']
    key = {col2[i]:i for i in range(0,len(col2))}
    for line in text:
        out = line.split('\t')
        # parses for other info
        info = [0]*len(key)
        for i in range(len(col1),len(out)):
            x = out[i].split(':')
            if x[0] in key:
                info[key[x[0]]] = x[-1]
        # save the data
        data.append(out[:len(col1)]+info)
    data = pd.DataFrame(data, columns=col1+col2)
    data = data.rename(columns={'cg':'CIGAR'})
    return data

def align_reads():
    '''
    This function aligns merged reads against the center sequences of each cluster using minimap2
    '''
    # run the alignment
    files = glob.glob('merged/*.fq*')
    files = [i.split('_')[0] for i in files]
    read1 = files
    
    # define results folder
    ofile = 'tmp/results.paf'
    df = []
    for i in range(len(read1)):
        print('aligning',read1[i], i, len(read1))
        c = 'minimap2 -c -x map-ont data/cluster_lin.fa '+read1[i]+' -o '+ofile
        subprocess.call(c, shell=True)
        x = parse_PAF(ofile)
        x['sample'] = read1[i].split('/')[-1].split('.')[0]
        
        # compute similarity score
        x['match'] = x['match'].astype(int)
        x['match_score'] = x['match'].astype(int)/x['t_len'].astype(int)
        
        # sort to best match
        x = get_best(x,'query_id')
        df.append(x)
    df = pd.concat(df)
    return df

def get_best(df, col, metric='match', stat='idxmax'):
    '''
    This function assigns each query to its best match
    '''
    df=df.reset_index(drop=True)
    idx = df.groupby(by=col).agg({metric:stat}).reset_index()
    return df.iloc[idx[metric].values].reset_index()

def get_feature_matrix(df):
    # get counts for each sample and feature
    col = ['sample','database_id','count']
    df['count'] = 1
    x = df[col].groupby(['sample','database_id']).sum().reset_index()
    # sort features by most frequently occuring samples
    y = x[col].groupby(['database_id']).sum().reset_index()
    y = y.sort_values(by='count')
    db = y.iloc[::-1]['database_id'].values
    jdict = {db[i]:i for i in range(len(db))}    
    # populate the indices for sample
    sample = x['sample'].drop_duplicates().values
    # initialize the feature matrix
    L = len(sample)
    N = len(db)
    mat = np.zeros([L,N])
    for i in range(len(sample)):
        x1 = x[x['sample']==sample[i]]
        for j,k in x1[['database_id','count']].values:
            j = jdict[j]
            mat[i,j] = k
        # normalize counts by total counts
        mat[i,:] = mat[i,:]/np.sum(mat[i,:])
    return mat
        
def get_taxon():
    # flag taxons
    # download silva Db to align 
    cmd = 'mmseqs databases SILVA databases/silvadb tmp'
    #subprocess.run(cmd.split(' '))
    
    # create query db
    cmd = 'mmseqs createdb clusters_rep_seq.fasta databases/queryDB'
    subprocess.run(cmd.split(' '))
    
    # search clusters against SILVA database
    cmd = 'mmseqs taxonomy databases/queryDB databases/silvadb result tmp'
    subprocess.run(cmd.split(' '))


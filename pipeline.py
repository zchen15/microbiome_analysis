import plasmid
import importlib
import glob
import subprocess

def merge_reads():
    # Collect files
    files = glob.glob('filtered/*')
    files = [i.split('_')[0] for i in files]
    # merge the reads
    for f in files:
        print('processing',f)
        r1 = f+'_1_filtered.fq.gz'
        r2 = f+'_2_filtered.fq.gz'
        fname = f.split('/')[-1]
        fname = 'merged/'+fname+'.fq.gz'
        aligner = plasmid.Aligner()
        aligner.merge_paired_end(r1, r2, fname)


def cluster_reads():
    # Cluster reads with linclust
    files = glob.glob('merged/*')
    cmd = 'mmseqs easy-linclust merged/*.fq.gz clusters tmp'
    subprocess.run(cmd.split(' '))

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

#merge_reads()
#cluster_reads()
get_taxon()
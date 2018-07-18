import sys
import pickle
import pprint 
"""
Script to generate concensus reads from R1 and R2 reads of deep sequencing data.

STEPs:
(BLAST data)
0. generate fastq
1. Open blast.out file and corresponding  .fastq files as Dictionary. {ID} : {...}
2. Extract IDs which have x.Frd in R1, and x.Rvs in R2.
3. For each of the extracted read, extract sequences.
4. Estimate the overlapping position by sliding the R2 sequence on R1, and define the best matched position.
5. Reconstitute the sequence based on qscore fir each position.
6. Generate new fasta file with R1 altered to concensus.
"""


def main(seq_dir,R1_id,R2_id):

    #R1_blast    = "%s/blast/out.blast/%s.blast"%(seq_dir,R1_id)
    #R2_blast    = "%s/blast/out.blast/%s.blast"%(seq_dir,R2_id)
    R1_fastq    = "%s/fragmented_fastq/%s.fastq"%(seq_dir,R1_id)
    R2_fastq    = "%s/fragmented_fastq/%s.fastq"%(seq_dir,R2_id)

    #r1_b = blast2dict(R1_blast)
    #r2_b = blast2dict(R2_blast)

    r1_f = fastq2dict(R1_fastq)
    r2_f = fastq2dict(R2_fastq)

    #ks   = r1_b.keys()
    #print len(r1_f.keys()), len(r2_f.keys())


    #Exrtracting only the reads that have target amplicons.
    c = 0
    filter_r1_f = r1_f
    filter_r2_f = r2_f
     
    out = []
    out2 = []
    ks =  filter_r1_f.keys() 
    #For each read, generating concensus read.
    for i in ks:
        r1_seq = filter_r1_f[i]["seq"]
        r1_q   = filter_r1_f[i]["qscore"]
        r2_seq = rev_comp(filter_r2_f[i]["seq"])
        r2_q   = filter_r2_f[i]["qscore"][::-1]
        #Sliding R1 over R2 to see the position where they overlap
        max_score = 0
        for j in range(5,len(r2_seq)-5):
            one = r1_seq[-j:]
            two = r2_seq[:j]
            score = seq_match(one,two)
            if score > max_score:
                pos       = j
                max_score = score 
                
        #print max_score, pos
        #print r2_seq
        #print "     ",r2_seq[:pos+5]
        add_seq = r1_seq[:-pos]
        #print add_seq
        for p in range(pos):
            if (qscore(r1_q[-pos+p]) >= qscore(r2_q[p]) ) :
                add_seq += r1_seq[-pos+p]
            else:
                add_seq += r2_seq[p]
            #print add_seq
            #print r1_q[-pos+p], r2_q[p], r1_seq[-pos+p],r2_seq[p]
        for p in range(len(r2_seq)-pos):
            add_seq += r2_seq[pos+p]
        l =  [i,add_seq]
        l2 = [i,filter_r2_f[i]["seq"]]
        out.append(l)
        out2.append(l2)
    #print out
    LL2fna(out,"%s/fragmented_fasta/%s.fna"%(seq_dir,R1_id))
    LL2fna(out2,"%s/fragmented_fasta/%s.fna"%(seq_dir,R2_id))
    #pickle.dump(reads_num,open("%s/QC/out.concensus/%s.pickle"%(seq_dir,R1_id) ,"wb"))

def seq_match(seq1,seq2):
    score = 0
    for i in range(len(seq1)):
        if (seq1[i] == seq2[i]):
            score +=1
    return score


def removekey(d, key):
    r = dict(d)
    del r[key]
    return r

def LL2fna(LL,name):
    with open(name,"w") as F:
        for L in LL:
            F.write(">%s\n" %  ("\n").join( [str(i) for i in L]))
    F.close()


def comp(seq):
    complement_dict = {'A':'T','T':'A','G':'C','C':'G',"N":"N"}
    return "".join([complement_dict[base] for base in seq])
def rev_comp(seq):
    complement_dict = {'A':'T','T':'A','G':'C','C':'G',"N":"N"}
    return "".join([complement_dict[base] for base in reversed(seq)])

def qscore(q):
    ref = {
            "!" : 0,
            '"' : 1,
            "#" : 2,
            "$" : 3,
            "%" : 4,
            "&" : 5,
            "'" : 6,
            "(" : 7,
            ")" : 8,
            "*" : 9,
            "+" : 10,
            "," : 11,
            "-" : 12,
            "." : 13,
            "/" : 14,
            "0" : 15,
            "1" : 16,
            "2" : 17,
            "3" : 18,
            "4" : 19,
            "5" : 20,
            "6" : 21,
            "7" : 22,
            "8" : 23,
            "9" : 24,
            ":" : 25,
            ";" : 26,
            "<" : 27,
            "=" : 28,
            ">" : 29,
            "?" : 30,
            "@" : 31,
            "A" : 32,
            "B" : 33,
            "C" : 34,
            "D" : 35,
            "E" : 36,
            "F" : 37,
            "G" : 38,
            "H" : 39,
            "I" : 40,
            "J" : 41,
            "K" : 42,
            "L" : 43,
            "M" : 44,
            "N" : 45,
            "O" : 46,
            "P" : 47,
            "Q" : 48,
            "R" : 49,
            "S" : 50,
            "T" : 51,
            "U" : 52,
            "V" : 53,
            "W" : 54,
            "X" : 55,
            "Y" : 56,
            "Z" : 57,
            "[" : 58}
    return ref[q]

def csv2LL(csv):
    LL = []
    with open(csv,"r") as F:
        for line in F:
            cols = line.split("\n")[0].split(",")
            LL.append(cols)
            #LL += cols
        F.close()
        return LL

def blast2dict(f_name):
    d = {}
    header = ["qseqid", "sseqid", "qstart", "qend", "sstart", "send", "evalue", "bitscore" "btop"]
    with open(f_name,"r") as F:
        c =0
        for line in F:
            c+=1
            cols = line.split(",")
            d[cols[0]] = {}
            for i in range(len(header[2:])+1):
                d[cols[0]][header[i+1].split("\n")[0]] =cols[i+1].split("\n")[0]

    return d

def fastq2dict(f_name):
    d = {}
    with open(f_name,"r") as F:
        c =0
        for line in F:
            c+=1
            if c == 1:
                nm  = line.split("\n")[0].split(">")[1]
            if c ==2:
                seq = line.split("\n")[0]
            if c ==3:
                qscore = line.split("\n")[0]
                c = 0
                d[nm] = {}
                d[nm]['seq'] = seq
                d[nm]['qscore'] = qscore

    return d


def load_data(pickle_or_multifasta):
    suffix =  pickle_or_multifasta.split(".")[-1]
    if (suffix == "pickle"):
        return pickle.load(open(pickle_or_multifasta))
    else:
        print "Input Error : The input is not a pickle file."


if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2],sys.argv[3])

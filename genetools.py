#When this script was written only god and I knew intimately how it works, now only god does
'''
This code analyses everything possible from a given DNA sequence, even if it's 'dirty' or pasted from a weird format. improvement is still possible, as always, but nothing beyond this point was created for readability so beware.
By all means use, reuse, distribute and hack away at this code.
I have no idea of the time complexity of this although I can testify that it runs the GFP sequence in a medium range smartphone without delay
'''
from textwrap import wrap
def translate(sequence):
    pseq=''
    for codon in sequence:
        if codon=='aug':
            pseq=pseq+'M'
        elif codon in ['gcu','gcc','gca','gcg']:
            pseq=pseq+'A'
        elif codon in ['cgu','cgc','cga','cgg','aga','agg']:
            pseq=pseq+'R'
        elif codon in ['aau','aac']:
            pseq=pseq+'N'
        elif codon in ['gau','gac']:
            pseq=pseq+'D'
        elif codon in ['ugu','ugc']:
            pseq=pseq+'C'
        elif codon in ['caa','cag']:
            pseq=pseq+'Q'
        elif codon in ['gaa','gag']:
            pseq=pseq+'E'
        elif codon in ['ggu','ggc','gga','ggg']:
            pseq=pseq+'G'
        elif codon in ['cau','cac']:
            pseq=pseq+'H'
        elif codon in ['auu','auc','aua']:
            pseq=pseq+'I'
        elif codon in ['cuu','cuc','cua','cug','uua','uug']:
            pseq=pseq+'L'
        elif codon in ['aaa','aag']:
            pseq=pseq+'K'
        elif codon in ['uuu','uuc']:
            pseq=pseq+'F'
        elif codon in ['ccu','ccc','cca','ccg']:
            pseq=pseq+'P'
        elif codon in ['ucu','ucc','uca','ucg','agu','agc']:
            pseq=pseq+'S'
        elif codon in ['acu','acc','aca','acge']:
            pseq=pseq+'T'
        elif codon=='ugg':
            pseq=pseq+'W'
        elif codon in ['uau','uac']:
            pseq=pseq+'Y'
        elif codon in ['guu','guc','gua','gug']:
            pseq=pseq+'V'
        elif codon in ['uaa','uga','uag']:
            break
    return pseq.upper()
def get_protein(seq):
    protseq=""
    rseq=""
    scod=None
    rseq=seq.replace('t','u')
    codseq=wrap(rseq,3)
    codseq1=wrap(rseq[1::],3)
    codseq2=wrap(rseq[2::],3)
    for cod in codseq:
        if cod=='aug':
            scod=cod
            print('the sequence is not shifted')
            shif=0
            break
        if scod==None:
            for cod in codseq1:
                if cod=='aug':
                    scod=cod
                    print('the sequence is shifted by 1')
                    shif=1
                    break
        if scod==None:
            for cod in codseq1:
                if cod=='aug':
                    scod=cod
                    print('the sequence is shifted by 2')
                    shif=2
                    break
    if scod!=None:
        print('this sequence codes for something!')
        if shif==0:
            protseq=translate(codseq)
        elif shif==1:
            protseq=translate(codseq1)
        elif shif==2:
            protseq=translate(codseq2)
    else:
            print('this sequence does not code for anything')
            protseq=None
    return tuple((protseq,rseq))
# uncomment this to try green fluorescent protein seq='1 atggtgagca agggcgagga gctgttcacc ggggtggtgc ccatcctggt cgagctggac       61ggcgacgtaa acggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctac    121ggcaagctgaccctgaagttcatctgcaccaccggcaagc tgcccgtgcc ctggcccacc181 ctcgtgacca ccctgaccta cggcgtgcag tgcttcagcc gctaccccga ccacatgaag  241 cagcacgact tcttcaagtc cgccatgccc gaaggctacg tccaggagcg caccatcttc301 ttcaaggacg acggcaacta caagacccgcgccgaggtga agttcgaggg cgacaccctg    361 gtgaaccgca tcgagctgaa gggcatcgacttcaaggagg acggcaacat cctggggcac     421 aagctagaat acaactacaacagccacaacgtctatatca tggccgacaagcagaagaac    481 ggcatcaagg tgaacttcaagatccgccac aacatcgagg acggcagcgtgcagctcgcc      541 gaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcc cgacaaccac      601tacctgagcacccagtccgc cctgagcaaagaccccaacgagaagcgcga tcacatggtc661 ctgctggagttcgtgaccgccgccgggatc actctcggcatggacgagct ctacaagtcc  721 tag'
seq=input('insert your sequence here: ').lower()
for char in seq:
    if char not in "acgt":
        seq=str(seq).replace(char, "")
ccode={'a':'t','g':'c','t':'a','c':'g'}
bcount={}
cseq=""
gcrat=None
palindrome=None
protein=''
rna=''
protweight=0
acon=None
acount={}
afrac={}
aweight={'G':57.021464,'P':97.052764,'A':71.037114,'V':99.068414,'L':113.084064,'I':113.084064,'M':131.040485,'C':103.009185,'F':147.068414,'Y':163.06332,'W':186.079313,'H':137.058912,'K':128.094963,'R':156.101111,'Q':128.058578,'N':114.042927,'E':129.042593,'D':115.026943,'S':87.032028,'T':101.047679}
for b in seq[::-1]:
    cseq=cseq+ccode[b]
for x in 'acgt':
    bcount[x]=seq.count(x)
gcrat=float('{:.3f}'.format((bcount['g']+bcount['c'])*100/len(seq)))
print('original sequence')
(protein,rna)=get_protein(seq)
if protein==None:
    print('inverted sequence')
    (protein,rna)=get_protein(seq[::-1])
if protein==None:
    print('complementary sequence')
    (protein,rna)=get_protein(cseq)
if protein==None:
    print('inverted complementary sequence')
    (protein,rna)=get_protein(cseq[::-1])
if protein!=None:
    for aa in protein:
        protweight+=aweight[aa]
    protweight=float('{:.4f}'.format(protweight+18.010564683))
    #TODO calculate the content and fraction of aminoacids
    for x in ''.join(set(protein)):
        acount[x]=protein.count(x)
        afrac[x]='{:.3}%'.format(100*acount[x]/len(protein))
    print('the sequence {} codes for {} is {} aminoacids long contains {} and weights {} daltons'.format(rna,protein,len(protein),set(protein),protweight))
    print(acount)
    print(afrac)
print(cseq)
print(bcount)
print('{}%'.format(gcrat))
if seq==cseq:
    print('{} is palindromic and its length is {}'.format(seq,len(seq)))
    palindrome=True
else:
    print('{} is not palindromic and its length is {}'.format(seq,len(seq)))
    palindrome=False

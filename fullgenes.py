def analyser(seq):
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
    print('{}%'.format(gcrat))
    if seq==cseq:
        print('{} is palindromic and its length is {}'.format(seq,len(seq)))
        palindrome=True
    else:
        print('{} is not palindromic and its length is {}'.format(seq,len(seq)))
        palindrome=False
    return tuple((palindrome,bcount,gcrat,cseq,rna,protein,protweight,acount,afrac))
seq1=str(input('insert the first sequence here: ')).lower()
for char in seq1:
    if char not in "acgt":
        seq1=str(seq1).replace(char, "")
seq2=str(input('insert the second sequence here: ')).lower()
for char in seq2:
    if char not in "acgt":
        seq2=str(seq2).replace(char, "")
if len(seq1)!=len(seq2):
    print('hamming distance can not be calculated')
    0/0
distance=0
charnum=0
transitions=0
transversions=0
purines='ag'
pirimidines='ct'
ttratio=None
for char in seq1:
    if char!=seq2[charnum]:
        distance+=1
        if ((char in purines) and (seq2[charnum] in pirimidines)) or ((char in pirimidines) and (seq2[charnum]) in purines):
            transversions+=1
        else:
            transitions+=1
    charnum+=1
if transversions!=0:
    ttratio=float('{:.3f}'.format(transitions/transversions))
else:
    ttratio='impossible to calculate'
print(transitions)
print(transversions)
print('the sequences have a hamming distance of {} and a ratio of transitions to transversions of {}'.format(distance,ttratio))
list(map(analyser,[seq1,seq2]))
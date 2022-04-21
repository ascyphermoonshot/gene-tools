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
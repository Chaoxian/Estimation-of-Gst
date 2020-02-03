import sys
from Bio.Seq import Seq
from Bio import SeqIO
from math import pow
from statistics import stdev
from statistics import mean

f1 = sys.argv[1]
f2 = open(sys.argv[2],'w')
S1 = int(sys.argv[3])
S2 = int(sys.argv[4])

All_list = []
All_gst, All_gst1, All_d = [],[],[]

for fasta in SeqIO.parse(f1, 'fasta'):
    All_list.append(fasta.seq)

Pop_s = len(All_list)
Seq_l = len(All_list[0])

# judge population size is right or not
if Pop_s != S1 + S2:
    print 'Wrong Population Size!!!'
    break

for i in range(Seq_1):
    All_site = []
    for j in range(Pop_s):
        All_site.append(All_list[j][i])
    site_spe = list(set((All_site)))
    Pop1 = All_site[:S1]
    Pop2 = All_site[S1:]
    Hs, Ht = 0.0,0.0
    if len(site_spe) == 1:
        continue
    elif len(site_spe) == 2:
        Alle1 = site_spe[0]
        Alle2 = site_spe[1]
        A1P1 = Pop1.count(Alle1)/float(S1)
        A1P2 = Pop2.count(Alle1)/float(S2)
        A2P1 = Pop1.count(Alle2)/float(S1)
        A2P2 = Pop2.count(Alle2)/float(S2)
        Hs = ((1 - pow(A1P1,2) - pow(A2P1,2)) + (1 - pow(A1P2,2) - pow(A2P2,2)))/2
        Ht = (1 - pow((A1P1 + A1P2)/2,2) - pow((A2P1 + A2P2)/2,2))
    elif len(site_spe) == 3:
        Alle1 = site_spe[0]
        Alle2 = site_spe[1]
        Alle3 = site_spe[2]
        A1P1 = Pop1.count(Alle1)/float(S1)
        A1P2 = Pop2.count(Alle1)/float(S2)
        A2P1 = Pop1.count(Alle2)/float(S1)
        A2P2 = Pop2.count(Alle2)/float(S2)
        A3P1 = Pop1.count(Alle3)/float(S1)
        A3P2 = Pop2.count(Alle3)/float(S2)
        Hs = ((1 - pow(A1P1,2) - pow(A2P1,2) - pow(A3P1,2)) + (1 - pow(A1P2,2) - pow(A2P2,2) - pow(A3P2,2)))/2
        Ht = (1 - pow((A1P1 + A1P2)/2,2) - pow((A2P1 + A2P2)/2,2) - pow((A3P1 + A3P2)/2,2))
    elif len(site_spe) == 4:
        Alle1 = site_spe[0]
        Alle2 = site_spe[1]
        Alle3 = site_spe[2]
        Alle4 = site_spe[3]
        A1P1 = Pop1.count(Alle1)/float(S1)
        A1P2 = Pop2.count(Alle1)/float(S2)
        A2P1 = Pop1.count(Alle2)/float(S1)
        A2P2 = Pop2.count(Alle2)/float(S2)
        A3P1 = Pop1.count(Alle3)/float(S1)
        A3P2 = Pop2.count(Alle3)/float(S2)
        A4P1 = Pop1.count(Alle4)/float(S1)
        A4P2 = Pop2.count(Alle4)/float(S2)
        Hs = ((1 - pow(A1P1,2) - pow(A2P1,2) - pow(A3P1,2) - pow(A4P1,2)) + (1 - pow(A1P2,2) - pow(A2P2,2) - pow(A3P2,2)) - pow(A4P2,2))/2
        Ht = (1 - pow((A1P1 + A1P2)/2,2) - pow((A2P1 + A2P2)/2,2) - pow((A3P1 + A3P2)/2,2) - pow((A4P1 + A4P2)/2))
    Gst_t = 1 - Hs/Ht
    Gst1_t = Gst_t/((1-Hs)/(1+Hs))
    D_t = 2*((Ht-Hs)/(1-Hs))
    All_gst.append(Gst_t)
    All_gst1.append(Gst1_t)
    All_d.append(D_t) 

Gst_m = mean(All_gst)
Gst_sd = stdev(All_gst)
Gst1_m = mean(All_gst1)
Gst1_sd = stdev(All_gst1)
D_m = mean(All_d)
D_sd = stdev(All_d)

print >> f2, 'Gst_m' + '\t' + 'Gst_sd' + '\t' + 'Gst1_m' + '\t' + 'Gst1_sd' + '\t' + 'D_m' + '\t' + 'D_sd'
print >> f2, str(Gst_m) + '\t' + str(Gst_sd) + '\t' + str(Gst1_m) + '\t' + str(Gst1_sd) + '\t' + str(D_m) + '\t' + str(D_sd)
f2.close()

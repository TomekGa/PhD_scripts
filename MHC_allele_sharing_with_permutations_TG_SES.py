#!/bin/python
#!@wiesiek 2018
#calculates alleles sharing and tests its significance
#modified by Tomek Gaczorek

#Input is tab-delimited text
#First line contains allele names (they should start in 3rd column), 1st and 2nd column have to have headers (e.g., ID Info)
#Other lines contain: c1 - individualID(not used by the script, so can be replaced with tab),
#c2 - Info about population: fields seperated with "_": species_zone_type_locality
#zone "0" is reserved for a situation when the same allopatric populations are used with all zone - only "0_allo" is used by ReadInAlleles
#type is "allo", "para", "syn'
#locality cannot include white spaces
#other fields may be present in c2, but are not used
#c3 and ff - presence absence of alleles encoded as 1/0

# running
#python3 <script> <input> <number_of_permutations> 

import random
import sys
import statistics

#for each zone, species, type reads individuals form a list of sets
#each set conatins alleles of a single individual
def ReadIndAlleles(In):
    Lcount=0
    #key: zone, value: dict of species within the zone;
    #zone "0" is reserved for allopatric pops common to multiple zones
    Zones={}
    InFile=open(In, 'r')
    for Line in InFile:
        Lcount += 1
        #print(Lcount)
        if Lcount == 1:
            SplitLine = Line.split()
            #loads allele names into list
            All=SplitLine[2:]
        #iterates through individuals and populates zones, species, and types
        else:
            SplitLine = Line.split()
            #gets individual's allele presence/absence into list
            Gen = SplitLine[2:]
            #reads zone, species, type
            Info = SplitLine[1].split("_")
            zone = Info[1]
            sp = Info[0]
            ty = Info[2]
            #initializes empty set for individual's alleles
            s = set()
            #populates the set with allele names
            i = 0
            for g in Gen:
                if g == '1': s.add(All[i])
                i += 1
            #populates the nested dict
            if zone in Zones:
                if sp in Zones[zone]:
                    if ty in Zones[zone][sp]:
                        Zones[zone][sp][ty].append(s)
                    else:
                        Zones[zone][sp][ty]=[]
                        Zones[zone][sp][ty].append(s)
                else:
                    Zones[zone][sp]={}
                    Zones[zone][sp][ty]=[]
                    Zones[zone][sp][ty].append(s)
            else:
                Zones[zone]={}
                Zones[zone][sp]={}
                Zones[zone][sp][ty]=[]
                Zones[zone][sp][ty].append(s)                                
    InFile.close()
    #if zone "0" exists then populate each zone "allo" with species info for that zone
    if "0" in Zones:
        for zone in Zones:
            if zone!="0":
                for sp in Zones[zone]:
                    if sp in Zones["0"]:
                        Zones[zone][sp]["allo"] = Zones["0"][sp]["allo"]
    return Zones

#takes four lists of sets, return tuple with fractions and numbers of shared alleles and total number of alleles
def Shared(allo1, allo2, para1, para2):
    a1 = set().union(*allo1)
    a2 = set().union(*allo2)
    p1 = set().union(*para1)
    p2 = set().union(*para2)
    Ns = [len(i) for i in (para1, allo1, para2, allo2)]
    NSpAl = [len(i) for i in (p1, a1, p2, a2)]
    nshAllo = float(len(a1.intersection(a2)))
    nshPara = float(len(p1.intersection(p2)))
    ntotAllo = float(len(a1.union(a2)))
    ntotPara = float(len(p1.union(p2)))
    return (nshAllo/ntotAllo, nshPara/ntotPara, nshAllo, nshPara, ntotAllo, ntotPara, nshPara/NSpAl[0], nshAllo/NSpAl[1], nshPara/NSpAl[2], nshAllo/NSpAl[3], Ns[0], Ns[1], Ns[2], Ns[3])

#permute individuals between types within species
def Perm(allo1, allo2, para1, para2, nperm):
    per1 = allo1 + para1
    per2 = allo2 + para2
    null_d = []
    for i in range(nperm):
        for j in (per1, per2): random.shuffle(j)
        rand_a1 = per1[:len(allo1)]
        rand_p1 = per1[len(allo1):]
        rand_a2 = per2[:len(allo2)]
        rand_p2 = per2[len(allo2):]
        s = Shared(rand_a1, rand_a2, rand_p1, rand_p2)
        dif_sh = s[1] - s[0]
        null_d.append(dif_sh)
    return null_d

#one sided test - proportion of nulls with difference larger or equal than observed
def Pval(to_test, null_d):
    one_sided_p = float(sum(i >= to_test for i in null_d))/len(null_d)
    return one_sided_p


#Calculates and writes to output allele sharing and other stats
def Body(Zones, Out, nperm):
    OutFile=open(Out, 'w')
    OutFile.write('Zone\tSp-Sp2\ttype\tNAllSharedClose\tNAllSharedFar\tPercAllSharedClose\tPercAllSharedFar\tpVal\tNperm\tPercAllSharedCloseSp1\tPercAllSharedFarSp1\t\
                    PercAllSharedCloseSp2\tPercAllSharedFarSp2\tNCloseSp1\tNFarSp1\tNCloseSp2\tNFarSp2\tSES\n')
    for zone in sorted(Zones.keys()):
        if zone!="0":
            spec = sorted(Zones[zone].keys())
            n_sp = len(Zones[zone])
            for i in range(n_sp):
                for j in range(n_sp):
                    if i < j:
                        sp1 = spec[i]
                        sp2 = spec[j]
                        #print (sp1, sp2)
                        #always test in order allo-para-sym
                        typ = sorted(Zones[zone][sp1].keys())
                        n_ty = len(Zones[zone][sp1])
                        for i in range(n_ty):
                            for j in range(n_ty):
                                if i < j:
                                    ty1 = typ[i]
                                    ty2 = typ[j]
                                    #print(ty1, ty2)
                                    a1 = Zones[zone][sp1][ty1]
                                    a2 = Zones[zone][sp2][ty1]
                                    p1 = Zones[zone][sp1][ty2]
                                    p2 = Zones[zone][sp2][ty2]
                                    obs = Shared(a1, a2, p1, p2)
                                    #print(a1)
                                    #print(a2)
                                    #print(p1)
                                    #print(p2)
                                    n_sh_ty1 = obs[2]
                                    n_sh_ty2 = obs[3]
                                    sh_ty1 = obs[0]
                                    sh_ty2 = obs[1]
                                    sh_ty1_sp1 = obs[6]
                                    sh_ty2_sp1 = obs[7]
                                    sh_ty1_sp2 = obs[8]
                                    sh_ty2_sp2 = obs[9]
                                    N_ty1_sp1 = obs[10]
                                    N_ty2_sp1 = obs[11]
                                    N_ty1_sp2 = obs[12]
                                    N_ty2_sp2 = obs[13]
                                    #print(sh_ty1)
                                    obs_diff = sh_ty2 - sh_ty1
                                    #print(obs_diff)
                                    null_d = Perm(a1, a2, p1, p2, nperm)
                                    #for i in range(len(null_d)): print(null_d[i])
                                    p = Pval(obs_diff, null_d)
                                    sd = statistics.stdev(null_d)
                                    if sd == 0:
                                        sd = 0.0001
                                    SES = (obs_diff - statistics.mean(null_d))/ sd
                                    OutFile.write('%s\t%s-%s\t%s-%s\t%d\t%d\t%4.1f\t%4.1f\t%5.4f\t%d\t%4.1f\t%4.1f\t%4.1f\t%4.1f\t%d\t%d\t%d\t%d\t%.4f\n'
                                                  % (zone, sp1, sp2, ty1, ty2, n_sh_ty2, n_sh_ty1, 100*sh_ty2, 100*sh_ty1, p, nperm,
                                                     100*sh_ty1_sp1, 100*sh_ty2_sp1, 100*sh_ty1_sp2, 100*sh_ty2_sp2, N_ty1_sp1, N_ty2_sp1, N_ty1_sp2, N_ty2_sp2,SES))
                        #OutFile.write('%s\t%sAllo-%sAllo\t%d\t%4.1f\t%d\t%d\n' % (zone, sp1, sp2, nshAllo, 100*(float(nshAllo)/ntotAllo), len(sp1allo), len(sp2allo)))
    OutFile.close()


zones = ReadIndAlleles(sys.argv[1])
outer_splited = sys.argv[1].split(".tsv")
outer = "".join([outer_splited[0],"_output.tsv"])
Body(zones, outer, int(sys.argv[2]))


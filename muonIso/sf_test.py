#!/usr/bin/env python
import os, sys, re
from ROOT import *
from rootUtil import useAtlasStyle, waitRootCmd, savehistory

gROOT.LoadMacro("HistFitting.C+")
from ROOT import effFitter

funlist=[]

def useEntryList(ch, cut, eName='elist', useExist=True):
    elist = gDirectory.Get(eName)

    if elist and not useExist:
        gDirectory.DeleteAll(eName)
        elist = None

    if  not elist:
        ch.Draw('>>'+eName, cut, 'entrylist')
        elist = gDirectory.Get(eName)

#     elist.Print()

    ### try to fix https://sft.its.cern.ch/jira/browse/ROOT-2739
    next1 = TIter(elist.GetLists());
    l1 = next1()
    while l1:
        l1.SetTreeName(ch.GetName())
#         print '-->',l1.GetTreeName()
        l1 = next1()
    elistf = elist.Clone()
#     print ch.GetName()
#     print "after fix"
#     elistf.Print()

#     for l in elist.GetLists():
#         l.SetTreeName(ch.GetName())

    ch.SetEntryList(elistf)


def test():
    dir0 = '/net/s3_data_home/dzhang/links/SAMPLES/R20/tpNtuple/v28'

    print dir0
    ch1 = TChain("ZmumuTPMuon/Trees/MuonProbe_OC_LooseProbes_noProbeIP/TPTree_MuonProbe_OC_LooseProbes_noProbeIP")
    ch2 = TChain("ZmumuTPMuon/Trees/MuonProbe_SC_LooseProbes_noProbeIP/TPTree_MuonProbe_SC_LooseProbes_noProbeIP")
#     ch1 = TChain('treeOS')
#     ch2 = TChain('treeSS')
#     files = dir0+'/bo_a_data16_*_v028*.root'
#     files = dir0+'/bo_a_data16_periodC*_v028*.root/ZmumuTPMuon/Trees/MuonProbe_OC_LooseProbes_noProbeIP/TPTree_MuonProbe_OC_LooseProbes_noProbeIP'
    files = dir0+'/bo_a_data16_periodC*_v028*.root'
    files1 = dir0+'/bo_a_data16_periodB*_v028*.root'
    ch1.Add(files)
    ch2.Add(files)
    ch1.Add(files1)
    ch2.Add(files1)

    ch3 = TChain("ZmumuTPMuon/Trees/MuonProbe_OC_LooseProbes_noProbeIP/TPTree_MuonProbe_OC_LooseProbes_noProbeIP")
    ch3.Add(dir0+"/a_mc15_361107_e3601_r7725_r7676_v028.root")

#     cut1 = 'probe_pt>35000&&probe_pt<40000'
#     cut1 = 'probe_pt>15000&&probe_pt<20000'
    cut1 = '(tag_matched_HLT_mu24_ivarmedium==1||tag_matched_HLT_mu26_ivarmedium==1)'
    cut1 += '&&probe_quality < 2'
    cut1 += '&&tag_ptcone40<1&&dilep_dphi>2'
    cut1 += '&&probe_pt>5000&&probe_pt<7000'
    useEntryList(ch1, cut1, 'elist1_lt20')
    useEntryList(ch2, cut1, 'elist2_lt20')
    useEntryList(ch3, cut1+'&&probe_truth_type==6&&tag_truth_type', 'elist3_lt20')
    ch3.Show(0)
#     useEntryList(ch3, cut1+'&&probe_truth_type==6', 'elist3_lt20')
#     ch1.Draw('>>elist1', cut1, 'entrylist')
#     elist1 = gDirectory.Get('elist1')
#     if elist1: print "elist1 exit" 
#     elist1.Print('v')
#     ch1.SetEntryList(elist1)
# 
#     ch2.Draw('>>elist2', cut1, 'entrylist')
#     elist2 = gDirectory.Get('elist2')
#     ch2.SetEntryList(elist2)

#     ch1.Draw("probe_q*tag_q:probe_matched_IsoTight")
    ch1.Draw("probe_pt*0.001>>(100,0,100)")
    ch2.Draw("probe_pt*0.001","","same")
    ch3.SetLineColor(6)
#     ch3.Draw("probe_pt*0.001","","same")
    waitRootCmd()


    h1 = TH1F("h1","h1;m_{#mu#mu} [GeV]",100,70,110)
    h_OSP = h1.Clone("h_OSP")
    ch1.Draw("dilep_mll*0.001>>h_OSP","probe_q*tag_q<0&&probe_matched_IsoTight==1")

    h_OSF = h1.Clone("h_OSF")
    ch1.Draw("dilep_mll*0.001>>h_OSF","probe_q*tag_q<0&&probe_matched_IsoTight==0")

    h_SSP = h1.Clone("h_SSP")
    ch2.Draw("dilep_mll*0.001>>h_SSP","probe_q*tag_q>0&&probe_matched_IsoTight==1")

    h_SSF = h1.Clone("h_SSF")
    ch2.Draw("dilep_mll*0.001>>h_SSF","probe_q*tag_q>0&&probe_matched_IsoTight==0")

    # MC
    h_MCP = h1.Clone("h_MCP")
    ch1.Draw("dilep_mll*0.001>>h_MCP","probe_q*tag_q<0&&probe_matched_IsoTight==1")

    h_MCF = h1.Clone("h_MCF")
    ch1.Draw("dilep_mll*0.001>>h_MCF","probe_q*tag_q<0&&probe_matched_IsoTight==0")



#     h_OSF.Draw()
#     h_SSF.SetLineColor(2)
#     h_SSF.Draw("same")
#     waitRootCmd()

    ef1 = effFitter()
    ef1.setMC(h_MCP, h_MCF)
    ef1.getEff(h_OSP, h_SSP, h_OSF, h_SSF, "test1")
    ef1.showValues()

    ### need to get histograms

    



funlist.append(test)

if __name__ == '__main__':
    savehistory('.')
    useAtlasStyle()
    for fun in funlist: print fun()

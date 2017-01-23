#!/usr/bin/env python
import os, sys, re
from ROOT import *
from rootUtil import useAtlasStyle, waitRootCmd, savehistory

gROOT.LoadMacro("HistFitting.C+")
from ROOT import effFitter

funlist=[]

def useEntryList(ch, cut, eName='elist', useExist=True, withfix=True):
    elist = gDirectory.Get(eName)

    if elist and not useExist:
        gDirectory.DeleteAll(eName)
        elist = None

    if  not elist:
        ch.Draw('>>'+eName, cut, 'entrylist')
        elist = gDirectory.Get(eName)

    ### try to fix https://sft.its.cern.ch/jira/browse/ROOT-2739
    if withfix:
        next1 = TIter(elist.GetLists());
        l1 = next1()
        while l1:
            l1.SetTreeName(ch.GetName())
            l1 = next1()
        elistf = elist.Clone()
        ch.SetEntryList(elistf)
        return
    else:
        ch.SetEntryList(elist)


class effChecker:
    def __init__(self, ch_data, ch_MC):
        self.ch1 = ch_data
        self.ch2 = None
        self.ch3 = ch_MC
        self.cacheFile = None
        self.cut0 = ''
        
    def updateEventList(self, cut1, elistTag=''):
        useEntryList(self.ch1, cut1, 'elist1'+elistTag)
        if self.ch2: useEntryList(self.ch2, cut1, 'elist2'+elistTag)
        useEntryList(self.ch3, cut1+'&&probe_truth_type==6&&tag_truth_type==6', 'elist3'+elistTag)

    def getHistosFromCache(self,saveKey):
        if self.cacheFile == None:
            return None

        h_OSP = self.cacheFile.Get("h_"+saveKey+"_OSP")
        h_OSF = self.cacheFile.Get("h_"+saveKey+"_OSF")
        h_SSP = self.cacheFile.Get("h_"+saveKey+"_SSP")
        h_SSF = self.cacheFile.Get("h_"+saveKey+"_SSF")
        h_MCP = self.cacheFile.Get("h_"+saveKey+"_MCP")
        h_MCF = self.cacheFile.Get("h_"+saveKey+"_MCF")

        return (h_MCP, h_MCF, h_OSP, h_OSF, h_SSP, h_SSF)
       
    def getHistosFromTree(self, ptRange, cache=True):
        elist0 = self.ch1.GetEventList()
        if elist0 == None and self.cut0 != '':
            self.updateEventList(self.cut0)

        saveKey = 'pt_%s_%s' % (str(ptRange[0]), str(ptRange[1]))
        cut1 = '&&probe_pt*0.001>%d&&probe_pt*0.001<%d' % ptRange

        h1 = TH1F("h1","h1;m_{#mu#mu} [GeV]",100,70,110)

        h_OSP = h1.Clone("h_OSP")
        self.ch1.Draw("dilep_mll*0.001>>h_OSP","probe_q*tag_q<0&&probe_matched_IsoTight==1"+cut1)

        h_OSF = h1.Clone("h_OSF")
        self.ch1.Draw("dilep_mll*0.001>>h_OSF","probe_q*tag_q<0&&probe_matched_IsoTight==0"+cut1)

        chx = self.ch2 if self.ch2!=None else self.ch1
        h_SSP = h1.Clone("h_SSP")
        chx.Draw("dilep_mll*0.001>>h_SSP","probe_q*tag_q>0&&probe_matched_IsoTight==1"+cut1)

        h_SSF = h1.Clone("h_SSF")
        chx.Draw("dilep_mll*0.001>>h_SSF","probe_q*tag_q>0&&probe_matched_IsoTight==0"+cut1)

        # MC
        h_MCP = h1.Clone("h_MCP")
        self.ch3.Draw("dilep_mll*0.001>>h_MCP","probe_q*tag_q<0&&probe_matched_IsoTight==1"+cut1)

        h_MCF = h1.Clone("h_MCF")
        self.ch3.Draw("dilep_mll*0.001>>h_MCF","probe_q*tag_q<0&&probe_matched_IsoTight==0"+cut1)

        if cache:
            self.cacheFile.cd()
            h_OSP.Write("h_"+saveKey+"_OSP")
            h_OSF.Write("h_"+saveKey+"_OSF")
            h_SSP.Write("h_"+saveKey+"_SSP")
            h_SSF.Write("h_"+saveKey+"_SSF")
            h_MCP.Write("h_"+saveKey+"_MCP")
            h_MCF.Write("h_"+saveKey+"_MCF")

        return (h_MCP, h_MCF, h_OSP, h_OSF, h_SSP, h_SSF)

    def test(self, ptRange = (7,9)):
        saveKey = 'pt_%s_%s' % (str(ptRange[0]), str(ptRange[1]))
        histos = self.getHistosFromCache(saveKey)
        if histos==None or histos[0]==None:
            histos = self.getHistosFromTree(ptRange)

        ef1 = effFitter()
        ef1.setMC(histos[0], histos[1])
        ef1.getEff(histos[2], histos[4], histos[3], histos[5], "test1")
        ef1.showValues()

    def test2(self, ptRange = (7,9)):
        saveKey = 'pt_%s_%s' % (str(ptRange[0]), str(ptRange[1]))
        histos = self.getHistosFromCache(saveKey)
        if histos==None or histos[0]==None:
            histos = self.getHistosFromTree(ptRange)

        ef1 = effFitter()
        ef1.setMC(histos[0], histos[1])
        ef1.hOS1 = histos[2]
        ef1.hSS1 = histos[3]
        ef1.hOS2 = histos[4]
        ef1.hSS2 = histos[5]
        ef1.eff = 0.9
        ef1.TF = 1
        ef1.showHists('test_x')
#         for i in range(10):
#             ef1.eff = 0.1*i
# #         ef1.TF = 1
#             ef1.showHists('test_'+str(i))



def test3():
    dir0 = '/net/s3_data_home/dzhang/links/SAMPLES/R20/tpNtuple/v28'

    print dir0
    ch1 = TChain("ZmumuTPMuon/Trees/MuonProbe_OC_LooseProbes_noProbeIP/TPTree_MuonProbe_OC_LooseProbes_noProbeIP")
    ch2 = TChain("ZmumuTPMuon/Trees/MuonProbe_SC_LooseProbes_noProbeIP/TPTree_MuonProbe_SC_LooseProbes_noProbeIP")
    files = dir0+'/bo_a_data16_periodC*_v028*.root'
    files1 = dir0+'/bo_a_data16_periodB*_v028*.root'
    ch1.Add(files)
    ch2.Add(files)
    ch1.Add(files1)
    ch2.Add(files1)

    ch3 = TChain("ZmumuTPMuon/Trees/MuonProbe_OC_LooseProbes_noProbeIP/TPTree_MuonProbe_OC_LooseProbes_noProbeIP")
    ch3.Add(dir0+"/a_mc15_361107_e3601_r7725_r7676_v028.root")

    cut1 = '(tag_matched_HLT_mu24_ivarmedium==1||tag_matched_HLT_mu26_ivarmedium==1)'
    cut1 += '&&probe_quality < 2'
    cut1 += '&&tag_ptcone40<1&&dilep_dphi>2'
    useEntryList(ch1, cut1, 'elist1_lt20')
    useEntryList(ch2, cut1, 'elist2_lt20')
#     useEntryList(ch3, cut1+'&&probe_truth_type==6&&tag_truth_type==6', 'elist3_lt20')
    useEntryList(ch3, cut1, 'elist3_lt20')

    ec1 = effChecker(ch1, ch3)
    ec1.cacheFile = TFile("testCache.root",'update')
    ec1.test()

# funlist.append(test3)


def test2():
    dir0 = '/net/s3_data_home/dzhang/work/muons/isolation/iso_20_7/lRun/output'

    print dir0
    ch1 = TChain("isoCutResults")
    ch1.Add(dir0+"/d16a1_Jan20a_tpOut.root")
    ch1.AddFriend("tree0")

    ch3 = TChain("isoCutResults")
    ch3.Add(dir0+"/m16a1_Jan20a_tpOut.root")
    ch3.AddFriend("tree0")

    cut1 = 'passGRL==1&&(tag_matched_HLT_mu24_ivarmedium==1||tag_matched_HLT_mu26_ivarmedium==1)'
    cut1 += '&&probe_quality < 2'
    cut1 += '&&fabs(tag_d0Sig)<3&&fabs(tag_z0C)<0.1&&fabs(probe_d0Sig)<3&&fabs(probe_z0C)<3'
    cut1 += '&&tag_ptcone40<1&&dilep_dphi>2'

#     cut1 += '&&probe_pt>5000&&probe_pt<7000'
#     useEntryList(ch1, cut1, 'elist1_lt20')
#     useEntryList(ch3, cut1+'&&probe_truth_type==6&&tag_truth_type==6', 'elist3_lt20')
#     ch3.Show(0)

#     ch1.Draw("probe_pt*0.001>>(100,0,100)","probe_q*tag_q<0","norm")
#     ch1.Draw("probe_pt*0.001","probe_q*tag_q>0","normsame")
#     ch3.SetLineColor(6)
#     ch3.Draw("probe_pt*0.001","probe_q*tag_q<0","normsame")
#     waitRootCmd()

    ec1 = effChecker(ch1, ch3)
    ec1.cut0 = cut1
    ec1.cacheFile = TFile("testCache1.root",'update')
    ec1.test()
#     ec1.test((30,40))

#     ec1.test2((30,40));

funlist.append(test2)


def test1():
    dir0 = '/net/s3_data_home/dzhang/work/muons/isolation/iso_20_7/lRun/output'

    print dir0
    ch1 = TChain("isoCutResults")
    ch1.Add(dir0+"/d16a1_Jan20a_tpOut.root")
    ch1.AddFriend("tree0")

#     tr0 = ch1.GetFriend("tree0")
#     tr0.SetBranchStatus("*truth*", 0)
#     tr0.SetBranchStatus("*Truth*", 0)

    ch3 = TChain("isoCutResults")
    ch3.Add(dir0+"/m16a1_Jan20a_tpOut.root")
    ch3.AddFriend("tree0")

    cut1 = 'passGRL==1&&(tag_matched_HLT_mu24_ivarmedium==1||tag_matched_HLT_mu26_ivarmedium==1)'
    cut1 += '&&probe_quality < 2'
    cut1 += '&&fabs(tag_d0Sig)<3&&fabs(tag_z0C)<0.1&&fabs(probe_d0Sig)<3&&fabs(probe_z0C)<3'
    cut1 += '&&tag_ptcone40<1&&dilep_dphi>2'

    cut1 += '&&probe_pt>5000&&probe_pt<7000'
    useEntryList(ch1, cut1, 'elist1_lt20')
    useEntryList(ch3, cut1+'&&probe_truth_type==6&&tag_truth_type==6', 'elist3_lt20')
    ch3.Show(0)

    ch1.Draw("probe_pt*0.001>>(100,0,100)","probe_q*tag_q<0","norm")
    ch1.Draw("probe_pt*0.001","probe_q*tag_q>0","normsame")
    ch3.SetLineColor(6)
    ch3.Draw("probe_pt*0.001","probe_q*tag_q<0","normsame")
    waitRootCmd()


    h1 = TH1F("h1","h1;m_{#mu#mu} [GeV]",100,70,110)
    h_OSP = h1.Clone("h_OSP")
    ch1.Draw("dilep_mll*0.001>>h_OSP","probe_q*tag_q<0&&probe_matched_IsoTight==1")

    h_OSF = h1.Clone("h_OSF")
    ch1.Draw("dilep_mll*0.001>>h_OSF","probe_q*tag_q<0&&probe_matched_IsoTight==0")

    h_SSP = h1.Clone("h_SSP")
    ch1.Draw("dilep_mll*0.001>>h_SSP","probe_q*tag_q>0&&probe_matched_IsoTight==1")

    h_SSF = h1.Clone("h_SSF")
    ch1.Draw("dilep_mll*0.001>>h_SSF","probe_q*tag_q>0&&probe_matched_IsoTight==0")

    # MC
    h_MCP = h1.Clone("h_MCP")
    ch3.Draw("dilep_mll*0.001>>h_MCP","probe_q*tag_q<0&&probe_matched_IsoTight==1")

    h_MCF = h1.Clone("h_MCF")
    ch3.Draw("dilep_mll*0.001>>h_MCF","probe_q*tag_q<0&&probe_matched_IsoTight==0")


#     h_OSF.Draw()
#     h_SSF.SetLineColor(2)
#     h_SSF.Draw("same")
#     waitRootCmd()

    ef1 = effFitter()
    ef1.setMC(h_MCP, h_MCF)
    ef1.getEff(h_OSP, h_SSP, h_OSF, h_SSF, "test1")
    ef1.showValues()

    ### need to get histograms
# funlist.append(test1)



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

    



# funlist.append(test)

if __name__ == '__main__':
    savehistory('.')
    useAtlasStyle()
    for fun in funlist: print fun()

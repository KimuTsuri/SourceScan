from ROOT import TH2F, gStyle, TCanvas
import ROOT
import json
import numpy as np
import sys


def main():
    args = sys.argv
    f = open(args[1], 'r')
    json_data = json.load(f)
    pix_list = ["Col", "Enable", "Hitbus", "InjEn", "TDAC"]
    col_list = np.linspace(264., 399., 136.)
    row_list = np.linspace(1., 191., 191.)
    #print(col_list)

    output = ROOT.TFile(args[1]+'enable_map.root','RECREATE')
    h2D = TH2F("h2D", "", 136, 264, 400, 192, 0, 192)
    gStyle.SetPalette(1)
#    ["PixelConfig"][{"Col", "Enable", "Hitbus", "InjEn", "TDAC"},{}]
    for col in col_list:
        #print(col)
        #print("Col:{} \t".format(col))#json_data["RD53A"]["PixelConfig"][col]))
        colnum = json_data["RD53A"]["PixelConfig"][int(col)]
        for row in row_list:
            if colnum["Enable"][int(row)] == 1:
                h2D.Fill(col, row)
    output.Write()
    can = ROOT.TCanvas()
    h2D.Draw("colz")
    
    h2D.SetXTitle("col")
    h2D.SetYTitle("row")
    h2D.SetStats(0)
    can.SaveAs(args[1]+'.png')


if __name__=='__main__':
    main()


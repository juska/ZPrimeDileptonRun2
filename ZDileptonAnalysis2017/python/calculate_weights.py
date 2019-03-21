#chad harrington 2019
#execute: python python/calculate_weights.py -l 25.6
import os
import argparse
from files import mcTups, ttbar_wgts

def main() :
  parser = argparse.ArgumentParser()
  parser.add_argument( "-l", "--lumi", type=float, action='store', required=True, help="Luminosity" )
  args = parser.parse_args()

  calc( args.lumi )

def calc( lumi ) :
  lumi = float(lumi)

  fname = "mc_weights.txt"
  file = open( "tmp_"+fname, "w+" )
  file.write( "dataset lumi xs lumi*xs events weight\n" )

  for mc in mcTups :
    file.write( "%s %.3f %.3f %.3f %i %f\n" % (mc.name, lumi, mc.xs, lumi*mc.xs, mc.nEvt, lumi*mc.xs/mc.nEvt) )

  for name, xs in ttbar_wgts :
    set = '_'.join(name.split('_')[:-1])
    set += '_total' if set != 'TTbar0-700' else ''
    evts = getEvts(set)
    if evts == None : continue
    file.write( "%s %.3f %.3f %.3f %i %f\n" % (name, lumi, xs, lumi*xs, evts, lumi*xs/evts) )

  file.close()

  os.system( "cat tmp_" + fname +  " | column -t > " + fname )
  os.system( "rm -f tmp_" + fname )

def getEvts( name ) :
  for mc in mcTups :
    if name == mc.name : return mc.nEvt

if __name__ == "__main__":
    main()

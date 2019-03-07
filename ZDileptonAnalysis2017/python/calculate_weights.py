#chad harrington 2019
#execute: python python/calculate_weights.py -l 25.6
import os
import argparse
from files import mcTups

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

  file.close()

  os.system( "cat tmp_" + fname +  " | column -t > " + fname )
  os.system( "rm -f tmp_" + fname )

if __name__ == "__main__":
    main()

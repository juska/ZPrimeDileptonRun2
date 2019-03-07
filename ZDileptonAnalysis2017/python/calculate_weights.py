import os
import argparse

#name, cross-section, nEvents
datasets = ( ('TTjets',831760,153531390), ('dylow',18610000,39521230) )

def main() :

  parser = argparse.ArgumentParser()
  parser.add_argument( "-l", "--lumi", type=float, action='store', required=True, help="Luminosity" )
  args = parser.parse_args()

  fname = "mc_weights.txt"
  file = open( "tmp_"+fname, "w+" )
  file.write( "dataset lumi xs lumi*xs events weight\n" )

  for set in datasets :
    file.write( "%s %.3f %.3f %.3f %i %.3f\n" % (set[0], args.lumi, set[1], args.lumi*set[1], set[2], args.lumi*set[1]/set[2]) )

  file.close()

  os.system( "cat tmp_" + fname +  " | column -t > " + fname )
  os.system( "rm -f tmp_" + fname )

if __name__ == "__main__":
    main()

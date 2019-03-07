chad = 'root://131.225.204.161:1094//store/user/cihar29/2017Analysis_trees'
bahareh = 'root://131.225.204.161:1094//store/user/broozbah/ZPRIME_2017'

from collections import namedtuple

dataTuple = namedtuple('dataTuple', 'dir name lumi')

dataDict = {
  'mm':dataTuple( chad,    'muon',                  41.86 ),
  'em':dataTuple( chad,    'muon',                  41.86 ),
  'ee':dataTuple( bahareh, 'DoubleEG_2017_31Mar18', 41.53 )
           }

mcTuple = namedtuple('mcTuple', 'dir name xs nEvt')

mcTups = (
  mcTuple( chad,    'TTjets',                          831760,      153531390 ),
#  mcTuple( chad,    'TT_M700to1000',                   1,           39258853 ),
#  mcTuple( chad,    'TT_M1000toInf',                   1,           20684462 ),
  mcTuple( chad,    'TTsemilep',                       453570,      41221873  ),
  mcTuple( chad,    'wjets',                           61526700,    77700506  ),
  mcTuple( chad,    'ww',                              118700,      7791498   ),
  mcTuple( chad,    'wz',                              47100,       3928630   ),
  mcTuple( chad,    'zz',                              16500,       1949768   ),
  mcTuple( chad,    'qcd30to80_bcToE',                 405623400,   16030011  ),
  mcTuple( chad,    'qcd80to170_bcToE',                38104430,    15972037  ),
  mcTuple( chad,    'qcd170to250_bcToE',               2635813.32,  9831904   ),
  mcTuple( chad,    'qcd250toInf_bcToE',               711925.875,  9921358   ),
  mcTuple( bahareh, 'DYJetsToLL_M10to50',              18610000,    39521230  ),
#  mcTuple( bahareh, 'DYJetsToLL_M50',                  6225420,     96897977  ), xs different from last year
  mcTuple( bahareh, 'ST_t-channel_top_4f',             136020,      5865875   ),
  mcTuple( bahareh, 'ST_t-channel_antitop_4f',         80950,       3939990   ),
  mcTuple( bahareh, 'ST_tW_top_5f_inclusiveDecays_v2', 35850,       7794186   ),
  mcTuple( bahareh, 'ST_tW_antitop_5f_v2',             35850,       7363960   ),
  mcTuple( bahareh, 'QCD_Pt-50to80_MuEnrichedPt5',     437504100,   24068613  ),
  mcTuple( bahareh, 'QCD_Pt-80to120_MuEnrichedPt5',    106033664.8, 23248995  ),
  mcTuple( bahareh, 'QCD_Pt-120to170_MuEnrichedPt5',   25190515.14, 20774848  ),
  mcTuple( bahareh, 'QCD_Pt-300to470_MuEnrichedPt5',   797352.69,   17744779  ),
  mcTuple( bahareh, 'QCD_Pt-470to600_MuEnrichedPt5',   79025.53776, 24243589  ),
  mcTuple( bahareh, 'QCD_Pt-600to800_MuEnrichedPt5',   25095.05908, 16392140  ),
  mcTuple( bahareh, 'QCD_Pt-800to1000_MuEnrichedPt5',  4707.368272, 15694987  ),
  mcTuple( bahareh, 'QCD_Pt-1000toInf_MuEnrichedPt5',  1621.31692,  11596693  )
         )

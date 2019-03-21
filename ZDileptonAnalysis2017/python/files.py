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
  mcTuple( chad,    'TTbar0-700',                      731209,      1 ), #waiting on updated TTJets sample
  mcTuple( chad,    'TTbar700-1000_total',             80104.6,     1  ),
  mcTuple( chad,    'TTbar1000-inf_total',             20446.8,     1  ),
  mcTuple( chad,    'TTsemilep',                       453570,      41221873  ),
  mcTuple( chad,    'wjets',                           61526700,    77700506  ),
  mcTuple( chad,    'ww',                              118700,      7791498   ),
  mcTuple( chad,    'wz',                              47100,       3928630   ),
  mcTuple( chad,    'zz',                              16500,       1949768   ),
  mcTuple( bahareh, 'DYJetsToLL_M10to50',              18610000,    39521230  ),
  mcTuple( bahareh, 'DYJetsToLL_M50',                  6225420,     96897977  ),
  mcTuple( bahareh, 'ST_t-channel_top_4f',             136020,      5865875   ),
  mcTuple( bahareh, 'ST_t-channel_antitop_4f',         80950,       3939990   ),
  mcTuple( bahareh, 'ST_tW_top_5f_inclusiveDecays_v2', 35850,       7794186   ),
  mcTuple( bahareh, 'ST_tW_antitop_5f_v2',             35850,       7363960   ),
         )

'''
  mcTuple( chad,    'qcd30to80_bcToE',                 405623400,   16030011  ),
  mcTuple( chad,    'qcd80to170_bcToE',                38104430,    15972037  ),
  mcTuple( chad,    'qcd170to250_bcToE',               2635813.32,  9831904   ),
  mcTuple( chad,    'qcd250toInf_bcToE',               711925.875,  9921358   ),
  mcTuple( bahareh, 'QCD_Pt-50to80_MuEnrichedPt5',     437504100,   24068613  ),
  mcTuple( bahareh, 'QCD_Pt-80to120_MuEnrichedPt5',    106033664.8, 23248995  ),
  mcTuple( bahareh, 'QCD_Pt-120to170_MuEnrichedPt5',   25190515.14, 20774848  ),
  mcTuple( bahareh, 'QCD_Pt-300to470_MuEnrichedPt5',   797352.69,   17744779  ),
  mcTuple( bahareh, 'QCD_Pt-470to600_MuEnrichedPt5',   79025.53776, 24243589  ),
  mcTuple( bahareh, 'QCD_Pt-600to800_MuEnrichedPt5',   25095.05908, 16392140  ),
  mcTuple( bahareh, 'QCD_Pt-800to1000_MuEnrichedPt5',  4707.368272, 15694987  ),
  mcTuple( bahareh, 'QCD_Pt-1000toInf_MuEnrichedPt5',  1621.31692,  11596693  )
         )
'''

'''
TTbar0-700_topPtWeightUP  726395
TTbar700-1000_topPtWeightUP 83359.7
TTbar1000-inf_topPtWeightUP 22005
TTbar0-700_pdfUP  731083
TTbar700-1000_pdfUP 80176.6
TTbar1000-inf_pdfUP 20500.7
TTbar0-700_pdfDN  731350
TTbar700-1000_pdfDN 80022.6
TTbar1000-inf_pdfDN 20386.9
TTbar0-700_q2UP 727477
TTbar700-1000_q2UP  82453.8
TTbar1000-inf_q2UP  21829.3
TTbar0-700_q2DN 734688
TTbar700-1000_q2DN  77947.3
TTbar1000-inf_q2DN  19124.2
'''

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
  mcTuple( chad,    'TTbar0-700',                      731209,      134082521 ),
  mcTuple( chad,    'TTbar700-1000_total',             80104.6,     54645901  ),
  mcTuple( chad,    'TTbar1000-inf_total',             20446.8,     24746283  ),
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
  mcTuple( bahareh, 'QCD_Pt-1000toInf_MuEnrichedPt5',  1621.31692,  11596693  ),
  mcTuple( bahareh, 'ZprimeToTT_M500_W5_total',        1000,        154498    ),
  mcTuple( bahareh, 'ZprimeToTT_M500_W50_total',       1000,        203071    ),
  mcTuple( bahareh, 'ZprimeToTT_M500_W150_total',      1000,        203071    ),
  mcTuple( bahareh, 'ZprimeToTT_M750_W7p5_total',      1000,        184233    ),
  mcTuple( bahareh, 'ZprimeToTT_M750_W75_total',       1000,        196945    ),
  mcTuple( bahareh, 'ZprimeToTT_M750_W225_total',      1000,        204767    ),
  mcTuple( chad,    'ZprimeToTT_M1000_W10_total',      1000,        187722    ),
  mcTuple( chad_local,    'ZprimeToTT_M1000_W100_ext', 1000,        182192    ),
  mcTuple( chad_local,    'ZprimeToTT_M1000_W300_ext', 1000,        171353    ),
  mcTuple( chad_local,    'ZprimeToTT_M1250_W12p5_ext',1000,        164623    ),
  mcTuple( chad,    'ZprimeToTT_M1250_W125_total',     1000,        189350    ),
  mcTuple( chad,    'ZprimeToTT_M1250_W375_total',     1000,        191525    ),
  mcTuple( chad,    'ZprimeToTT_M1500_W15_total',      1000,        204773    ),
  mcTuple( chad,    'ZprimeToTT_M1500_W150_total',     1000,        191088    ),
  mcTuple( chad,    'ZprimeToTT_M1500_W450_total',     1000,        195127    ),
  mcTuple( chad,    'ZprimeToTT_M2000_W20_total',      1000,        210348    ),
  mcTuple( chad,    'ZprimeToTT_M2000_W200_total',     1000,        222006    ),
  mcTuple( chad,    'ZprimeToTT_M2000_W600_total',     1000,        194303    ),
  mcTuple( chad,    'ZprimeToTT_M2500_W25_total',      1000,        207452    ),
  mcTuple( chad,    'ZprimeToTT_M2500_W250_total',     1000,        199498    ),
  mcTuple( chad_local,    'ZprimeToTT_M2500_W750',     1000,        19748     ),
  mcTuple( bahareh, 'ZprimeToTT_M3000_W30_total',      1000,        189400    ),
  mcTuple( bahareh, 'ZprimeToTT_M3000_W300_total',     1000,        179744    ),
  mcTuple( bahareh, 'ZprimeToTT_M3000_W900_total',     1000,        213531    ),
  mcTuple( bahareh, 'ZprimeToTT_M3500_W35_total',      1000,        198848    ),
  mcTuple( bahareh, 'ZprimeToTT_M3500_W350_total',     1000,        209715    ),
  mcTuple( bahareh,    'ZprimeToTT_M3500_W1050_total', 1000,        196070    ),
  mcTuple( bahareh, 'ZprimeToTT_M4000_W40_total',      1000,         197552    ),
  mcTuple( bahareh,    'ZprimeToTT_M4000_W400_total',  1000,         206518    ),
  mcTuple( bahareh, 'ZprimeToTT_M4000_W1200_total',    1000,         200961    ),
  mcTuple( bahareh,    'ZprimeToTT_M5000_W50_total',    1000,       174358    ),
  mcTuple( bahareh, 'ZprimeToTT_M5000_W500_total.root', 1000,       206460    ),
  mcTuple( bahareh, 'ZprimeToTT_M5000_W1500_total.root',1000,       192519    ),
  mcTuple( bahareh, 'ZprimeToTT_M6000_W60_total',       1000,       200004     ),
  mcTuple( bahareh, 'ZprimeToTT_M6000_W600_total',      1000,       205469     ),
  mcTuple( bahareh, 'ZprimeToTT_M6000_W1800_total',     1000,       195251     ),
  mcTuple( bahareh,    'ZprimeToTT_M7000_W70_total',    1000,       171646     ),
  mcTuple( bahareh,    'ZprimeToTT_M7000_W700_total',   1000,       200125     ),
  mcTuple( bahareh,    'ZprimeToTT_M7000_W2100_total',  1000,       11713      ),
  mcTuple( bahareh,    'ZprimeToTT_M8000_W80_total',    1000,       187899     ),
  mcTuple( bahareh,    'ZprimeToTT_M8000_W800_total',   1000,       194374     ),
  mcTuple( bahareh,    'ZprimeToTT_M8000_W2400_total',  1000,       198477     ),
  mcTuple( bahareh,    'ZprimeToTT_M9000_W900_total',   1000,       220286     ),
  mcTuple( bahareh,    'ZprimeToTT_M9000_W2700_total',  1000,       191327     ),
  mcTuple( bahareh,    'RSGluonToTT_M-500',             1000,       200000     ),
  mcTuple( bahareh,    'RSGluonToTT_M-750',             1000,       184000     ),
  mcTuple( bahareh,    'RSGluonToTT_M-1000',            1000,       200000     ),
  mcTuple( bahareh,    'RSGluonToTT_M-1250',            1000,       192000     ),
  mcTuple( bahareh,    'RSGluonToTT_M-1500',            1000,       200000     ),
  mcTuple( bahareh,    'RSGluonToTT_M-2000',            1000,       200000     ),
  mcTuple( bahareh,    'RSGluonToTT_M-2500',            1000,       188000     ),
  mcTuple( bahareh,    'RSGluonToTT_M-3000',            1000,       200000     ),
  mcTuple( bahareh,    'RSGluonToTT_M-3500',            1000,       192000     ),
  mcTuple( bahareh,    'RSGluonToTT_M-4000',            1000,       184000     ),
  mcTuple( bahareh,    'RSGluonToTT_M-4500',            1000,       194000     ),
  mcTuple( bahareh,    'RSGluonToTT_M-5000',            1000,       197000     ),
  mcTuple( bahareh,    'RSGluonToTT_M-5500',            1000,       194000     ),
  mcTuple( bahareh,    'RSGluonToTT_M-6000',            1000,       200000     ),
  mcTuple( bahareh,    'RSGluonToTT_M-6500',            1000,       197000     ),
  mcTuple( bahareh,    'RSGluonToTT_M-7000',            1000,       196000     ),
  mcTuple( bahareh,    'RSGluonToTT_M-7500',            1000,       192000     ),
  mcTuple( bahareh,    'RSGluonToTT_M-8000',            1000,       200000     )
         )

ttbar_wgts = (
  ('TTbar0-700_topPtWeightUP',    726395  ),
  ('TTbar700-1000_topPtWeightUP', 83359.7 ),
  ('TTbar1000-inf_topPtWeightUP', 22005   ),
  ('TTbar0-700_pdfUP',            731083  ),
  ('TTbar700-1000_pdfUP',         80176.6 ),
  ('TTbar1000-inf_pdfUP',         20500.7 ),
  ('TTbar0-700_pdfDN',            731350  ),
  ('TTbar700-1000_pdfDN',         80022.6 ),
  ('TTbar1000-inf_pdfDN',         20386.9 ),
  ('TTbar0-700_q2UP',             727477  ),
  ('TTbar700-1000_q2UP',          82453.8 ),
  ('TTbar1000-inf_q2UP',          21829.3 ),
  ('TTbar0-700_q2DN',             734688  ),
  ('TTbar700-1000_q2DN',          77947.3 ),
  ('TTbar1000-inf_q2DN',          19124.2 )
             )
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
  mcTuple( chad,    'TTbar0-700',                      731209,      134082521 ),
  mcTuple( chad,    'TTbar700-1000_total',             80104.6,     54645901  ),
  mcTuple( chad,    'TTbar1000-inf_total',             20446.8,     24746283  ),
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

ttbar_wgts = (
  ('TTbar0-700_topPtWeightUP',    726395  ),
  ('TTbar700-1000_topPtWeightUP', 83359.7 ),
  ('TTbar1000-inf_topPtWeightUP', 22005   ),
  ('TTbar0-700_pdfUP',            731083  ),
  ('TTbar700-1000_pdfUP',         80176.6 ),
  ('TTbar1000-inf_pdfUP',         20500.7 ),
  ('TTbar0-700_pdfDN',            731350  ),
  ('TTbar700-1000_pdfDN',         80022.6 ),
  ('TTbar1000-inf_pdfDN',         20386.9 ),
  ('TTbar0-700_q2UP',             727477  ),
  ('TTbar700-1000_q2UP',          82453.8 ),
  ('TTbar1000-inf_q2UP',          21829.3 ),
  ('TTbar0-700_q2DN',             734688  ),
  ('TTbar700-1000_q2DN',          77947.3 ),
  ('TTbar1000-inf_q2DN',          19124.2 )
             )


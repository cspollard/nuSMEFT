import yoda
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

# neglects overflows
def numpy1d(h):
  return \
    np.array \
    (  [ h.bin(xidx).val()
          for xidx in range(1, h.numBinsX()+1)
        ]
      
    )


def integrate(arr):
  integral = 0
  for val in arr:
    integral += val
  return integral



# h_sm = yoda.read(argv[1])["/MyAnalysis/" + str(argv[3])]
# h_bsm = yoda.read(argv[2])["/MyAnalysis/" + str(argv[3])]
# edges = h_sm.xEdges()

# hist_sm = numpy1d(h_sm)
# hist_bsm = numpy1d(h_bsm)

# #normalise
# hist_sm = hist_sm / integrate(hist_sm)
# hist_bsm = hist_bsm / integrate(hist_bsm)

# plt.stairs(hist_sm, edges, label='sm')
# plt.stairs(hist_bsm, edges, label='bsm')
# plt.legend()
# plt.xticks(np.arange(min(edges),max(edges),1),minor=True)
# plt.title("Reconstructed $u_T$ before event selection (Atlas)")
# plt.savefig('/data/atlas/users/wadh6495/rivet/presentation/plots/eminus/' + str(argv[3]) + '.pdf')

plots = [("hist_e_pT","Reconstructed $e$ $p_T$"), \
          ("hist_mT", "Reconstructed $m_T$"), \
          ("hist_pT_miss", "Reconstructed missing $p_T$"), \
          ("hist_e_pT_true", "True $e$ $p_T$"), \
          ("hist_N_pT_true", "True neutrino $p_T$"), \
          ("hist_pT_miss_true", "True missing $p_T$"), \
          ("hist_mT_epT_miss_true", "Transverse mass calculated using true missing $p_T$ and true $e$ $p_T$"), \
          ("hist_mT_eN_true", "Transverse mass calculated using true neutrino $p_T$ and true $e$ $p_T$"), \
          ("hist_m_inv_true", "Invariant mass"), \
          ("hist_e_pT_before", "$e$ $p_T$ before event selection" ), \
          ("hist_mT_before", "Transverse mass before event selection"),\
          ("hist_pT_miss_before", "Missing $p_T$ before event selection"),\
          ("hist_e_abseta_before", "$| \eta |$ before event selection"),\
          ("hist_E_T_U","$u_T$ before event selection")]

for val in plots:
  h_sm = yoda.read(argv[1])["/MyAnalysis/" + val[0]]
  h_bsm = yoda.read(argv[2])["/MyAnalysis/" + val[0]]
  edges = h_sm.xEdges()

  hist_sm = numpy1d(h_sm)
  hist_bsm = numpy1d(h_bsm)

  #normalise
  hist_sm = hist_sm / integrate(hist_sm)
  hist_bsm = hist_bsm / integrate(hist_bsm)

  fig = plt.figure()

  plt.stairs(hist_sm, edges, label='sm')
  plt.stairs(hist_bsm, edges, label='bsm')
  plt.legend()
  plt.xticks(np.arange(min(edges),max(edges),1),minor=True)
  plt.title(val[1])
  plt.savefig(argv[3] + val[0] + '.pdf')
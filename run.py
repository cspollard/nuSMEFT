import yoda
import numpy
from matplotlib import figure
from sys import argv

scan = yoda.read(argv[1])["/MyAnalysis/hist_e_pT_2d"]

# neglects overflows
def numpy2d(h):
  return \
    numpy.array \
    ( [ [ h.bin(xidx, yidx).val()
          for yidx in range(1, h.numBinsY()+1)
        ]
        for xidx in range(1, h.numBinsX()+1)
      ]
    )

hists = numpy2d(scan)

edges = list(scan.yEdges())

fig = figure.Figure((6, 6))

plt = fig.add_subplot()

for i in range(hists.shape[0]):
  plt.stairs \
    ( hists[i]
    , edges = edges
    , fill = False
    )


fig.savefig("test.png")

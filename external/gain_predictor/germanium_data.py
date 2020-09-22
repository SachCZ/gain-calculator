"""calculate the electron impact excitation cross sections
"""
import sys
from pfac import fac


use_openmp = False
if len(sys.argv) == 2 and sys.argv[1] == 'openmp':
    use_openmp = True

if use_openmp:
    # enable openmp with 2 cores
    fac.InitializeMPI(2)

fac.SetAtom('Ge')
# 1s shell is closed
fac.Closed('1s 2s')
fac.Config('2p6', group = 'base')
fac.Config('2p5 3s1', group = 'lower')
fac.Config('2p5 3p1', group = 'upper')

# Self-consistent iteration for optimized central potential
fac.ConfigEnergy(0)
fac.OptimizeRadial(['base', 'lower', 'upper'])
fac.ConfigEnergy(1)
fac.Structure('ne.lev.b')
fac.MemENTable('ne.lev.b')
fac.PrintTable('ne.lev.b', 'ne.lev', 1)

fac.SetCEGrid(10, 10, 20000)
fac.CETable('ne.ce.b', ['base', 'lower'], ['lower', 'upper'])
fac.PrintTable('ne.ce.b', 'ne.ce', 1)

fac.TransitionTable('ne.tr.b', ['base', 'lower'], ['lower', 'upper'])
fac.PrintTable('ne.tr.b', 'ne.tr', 1)

if use_openmp:
    fac.FinalizeMPI()

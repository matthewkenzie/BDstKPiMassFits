# A module that helps convert known particle types into reasonable looking TLatex

Bp  = '#it{B}^{#plus}'
Bm  = '#it{B}^{#minus}'
Bz  = '#it{B}#lower[0.1]{^{#scale[0.85]{0}}}'
Bzb = '#bar{#it{B}}#kern[0.1]{}#lower[0.1]{^{#scale[0.85]{0}}}'
Bd  = Bz
Bdb = Bzb
Bs  = '#it{B}#lower[0.1]{^{#scale[0.85]{0}}}#kern[-1.4]{#lower[-0.2]{_{#it{s}}}}'
Bsb = '#bar{#it{B}}#lower[0.1]{^{#scale[0.85]{0}}}#kern[-1.4]{#lower[-0.2]{_{#it{s}}}}'
Lb  = '#it{#Lambda}#lower[0.1]{^{#scale[0.85]{0}}}#kern[-0.6]{#lower[-0.2]{_{#it{b}}}}'
Dz  = '#it{D}#kern[0.2]{}#lower[0.1]{^{#scale[0.85]{0}}}'
Dzb = '#bar{#it{D}}#kern[0.2]{}#lower[0.1]{^{#scale[0.85]{0}}}'
g   = '#it{#gamma}'
piz = '#it{#pi}#kern[0.1]{^{#scale[0.85]{0}}}'
pip = '#it{#pi}#kern[0.1]{^{#plus}}'
pim = '#it{#pi}#kern[0.1]{^{#minus}}'
Kz  = '#it{K}^{0}'
Kp  = '#it{K}^{#plus}'
Km  = '#it{K}^{#minus}'
p   = '#it{p}'

def Decay(m, d1=None, d2=None, d3=None, d4=None, d5=None):
  assert( m in globals().keys() )
  if d1: assert( d1 in globals().keys() )
  if d2: assert( d2 in globals().keys() )
  if d3: assert( d3 in globals().keys() )
  if d4: assert( d4 in globals().keys() )
  if d5: assert( d5 in globals().keys() )
  ret_str = globals()[m] + ' #rightarrow '
  if d1: ret_str += globals()[d1]
  if d2: ret_str += globals()[d2]
  if d3: ret_str += globals()[d3]
  if d4: ret_str += globals()[d4]
  if d5: ret_str += globals()[d5]
  return ret_str

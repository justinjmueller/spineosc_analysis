import uproot, numpy as np

WSYST = '/Users/mueller/data/spineosc/v1.1.0/spineosc_numu_sbnd_v1.1.0_wsyst.root'

f = uproot.open(WSYST)
t = f['events/sbnd_large/selected']
br = t.arrays(['reco_vertex_x', 'reco_leading_muon_end_x', 'reco_fiducialize_cathode',
               'reco_veto_sbnd_highy_highz', 'reco_leading_muon_muon_softmax'], library='np')

base = (br['reco_fiducialize_cathode'] == 1) & (br['reco_veto_sbnd_highy_highz'] == 1)
crosser = (br['reco_vertex_x'] * br['reco_leading_muon_end_x']) < 0

score = br['reco_leading_muon_muon_softmax']

sel_cross = base & crosser
sel_noncross = base & (~crosser)

n_cross = sel_cross.sum()
n_noncross = sel_noncross.sum()
print(f"baseline-selected: crossers={n_cross} ({100*n_cross/base.sum():.1f}%)  non-crossers={n_noncross} ({100*n_noncross/base.sum():.1f}%)")

edges = np.arange(0.0, 1.0001, 0.05)
h_cross, _ = np.histogram(score[sel_cross], bins=edges)
h_noncross, _ = np.histogram(score[sel_noncross], bins=edges)

h_cross_norm = h_cross / h_cross.sum()
h_noncross_norm = h_noncross / h_noncross.sum()

print(f"\nmean muon_score: crossers={np.mean(score[sel_cross]):.4f}  non-crossers={np.mean(score[sel_noncross]):.4f}")
print(f"median muon_score: crossers={np.median(score[sel_cross]):.4f}  non-crossers={np.median(score[sel_noncross]):.4f}")
print(f"frac with muon_score < 0.5: crossers={np.mean(score[sel_cross]<0.5):.4f}  non-crossers={np.mean(score[sel_noncross]<0.5):.4f}")
print(f"frac with muon_score < 0.8: crossers={np.mean(score[sel_cross]<0.8):.4f}  non-crossers={np.mean(score[sel_noncross]<0.8):.4f}")
print(f"frac with muon_score > 0.95: crossers={np.mean(score[sel_cross]>0.95):.4f}  non-crossers={np.mean(score[sel_noncross]>0.95):.4f}")

print(f"\n{'bin':>14}{'cross_n':>10}{'cross_frac':>12}{'noncross_n':>12}{'noncross_frac':>15}{'ratio(c/nc)':>13}")
for i in range(len(h_cross)):
    r = h_cross_norm[i]/h_noncross_norm[i] if h_noncross_norm[i] > 0 else float('nan')
    print(f"[{edges[i]:>5.2f},{edges[i+1]:>5.2f})  {h_cross[i]:>8d}  {h_cross_norm[i]:>10.4f}  {h_noncross[i]:>10d}  {h_noncross_norm[i]:>13.4f}  {r:>11.3f}")

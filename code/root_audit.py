import numpy as np,csv
from pathlib import Path
from scipy.linalg import eigvals
from lk_model import Params,solve_ecm,matrices,delta
from lk_analysis import characteristic_roots,chain_matrix
BASE=Path(__file__).resolve().parents[1]
(BASE/'new_results').mkdir(exist_ok=True)
(BASE/'data').mkdir(exist_ok=True)

def winding_count(p,e,left=-.05,right=.02,ymax=.5,n=1000):
    A0,A1=matrices(p,e)
    # positive orientation rectangle: bottom L->R, right bottom->top, top R->L, left top->bottom
    z=np.concatenate([
      np.linspace(left-1j*ymax,right-1j*ymax,n,endpoint=False),
      np.linspace(right-1j*ymax,right+1j*ymax,n,endpoint=False),
      np.linspace(right+1j*ymax,left+1j*ymax,n,endpoint=False),
      np.linspace(left+1j*ymax,left-1j*ymax,n+1)])
    vals=np.array([np.linalg.det(delta(x,A0,A1,p.tau)) for x in z])
    phase=np.unwrap(np.angle(vals))
    return int(np.rint((phase[-1]-phase[0])/(2*np.pi))),np.min(np.abs(vals))

def main():
 cases=[('solitary',Params(kappa=0,tau=20,phi=0)),('tradeoff',Params(kappa=.008,tau=20,phi=0)),('win_win',Params(kappa=.009,tau=20,phi=np.pi/4))]
 rows=[]
 for name,p in cases:
  e=solve_ecm(p);count,mind=winding_count(p,e)
  rr=characteristic_roots(p,e,M=100,nseeds=100,real_cut=-.2,imag_cut=.6)
  inbox=[z for z in rr if -.05<z.real<.02 and abs(z.imag)<.5]
  conv=[]
  for M in [30,45,60,90,120]:
   vals=eigvals(chain_matrix(p,e,M)); vals=vals[(np.abs(vals)>1e-7)&(vals.real>-.1)&(np.abs(vals.imag)<.5)]
   vals=vals[np.argsort(vals.real)[::-1]][:2]
   conv.append(vals)
  for M,vals in zip([30,45,60,90,120],conv):
   lead=vals[0] if len(vals) else np.nan+1j*np.nan
   rows.append([name,count,mind,len(inbox),M,lead.real,lead.imag])
  with open(BASE/'new_results'/f'roots_{name}.txt','w') as f:
   f.write(f'argument principle count={count}, min det on contour={mind:.3e}\n')
   f.write('refined roots in rectangle:\n'+'\n'.join(str(z) for z in inbox)+'\n')
 with open(BASE/'data'/'root_audit.csv','w',newline='') as f:
  w=csv.writer(f);w.writerow(['case','argument_count','min_det_contour','refined_roots_in_box','chain_M','chain_lead_real','chain_lead_imag']);w.writerows(rows)
 print(open(BASE/'data'/'root_audit.csv').read())
if __name__=='__main__':main()

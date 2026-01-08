# --- fa_volume_compare.py (RDKit  ensemble → geometry/volume summaries) ---
# Usage:  python fa_volume_compare.py --input "FA SMILES.csv" --out "FA_boltzmann_summaries.csv" --nconf 64 --grid 0.5
import argparse, json, sys, platform, datetime
import numpy as np, pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolTransforms import ComputePrincipalAxesAndMoments
PT = Chem.GetPeriodicTable()
def _rg(m,cid,include_h=True):
    conf=m.GetConformer(cid)
    idxs=range(m.GetNumAtoms()) if include_h else [i for i in range(m.GetNumAtoms()) if m.GetAtomWithIdx(i).GetAtomicNum()!=1]
    masses=np.array([PT.GetAtomicWeight(m.GetAtomWithIdx(i).GetAtomicNum()) for i in idxs],float)
    coords=np.array([conf.GetAtomPosition(i) for i in idxs],float)
    com=np.average(coords,axis=0,weights=masses)
    return float(np.sqrt(np.average(((coords-com)**2).sum(axis=1),weights=masses)))
def _extents(m,cid):
    conf=m.GetConformer(cid); axes,_=ComputePrincipalAxesAndMoments(conf,ignoreHs=True)
    A=np.array(axes); coords=np.array([conf.GetAtomPosition(i) for i in range(m.GetNumAtoms()) if m.GetAtomWithIdx(i).GetAtomicNum()!=1],float)
    proj=coords@A.T; L,W,T=(proj[:,k].max()-proj[:,k].min() for k in range(3)); return float(L),float(W),float(T)
def _end2end(m,cid):
    conf=m.GetConformer(cid); c0=None
    for a in m.GetAtoms():
        if a.GetAtomicNum()!=6: continue
        if any(n.GetAtomicNum()==8 and m.GetBondBetweenAtoms(a.GetIdx(),n.GetIdx()).GetBondTypeAsDouble()==2.0 for n in a.GetNeighbors()):
            c0=a.GetIdx(); break
    if c0 is None: c0=0
    p0=np.array(conf.GetAtomPosition(c0))
    terms=[a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum()==6 and a.GetDegree()==1 and not any(n.GetAtomicNum()==8 for n in a.GetNeighbors())]
    if not terms: terms=[a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum()==6]
    far=max(terms,key=lambda i: np.linalg.norm(np.array(conf.GetAtomPosition(i))-p0))
    return float(np.linalg.norm(np.array(conf.GetAtomPosition(far))-p0))
def _npr(m,cid):
    conf=m.GetConformer(cid); _,mom=ComputePrincipalAxesAndMoments(conf,ignoreHs=True); I=sorted(list(mom)); I1,I2,I3=I
    return float(I1),float(I2),float(I3), (float(I1/I3) if I3 else 0.0), (float(I2/I3) if I3 else 0.0)
def _grid_union(coords,radii,grid=0.5,margin=1.0):
    coords=np.asarray(coords,float); radii=np.asarray(radii,float)
    mins=coords.min(axis=0)-(radii.max()+margin); maxs=coords.max(axis=0)+(radii.max()+margin)
    xs=np.arange(mins[0],maxs[0]+grid*0.5,grid); ys=np.arange(mins[1],maxs[1]+grid*0.5,grid); zs=np.arange(mins[2],maxs[2]+grid*0.5,grid)
    vol=grid**3; nx,ny=len(xs),len(ys); occ=0
    for z in zs:
        plane=np.zeros((nx,ny),dtype=bool)
        for (x0,y0,z0),r in zip(coords,radii):
            dz2=(z-z0)**2
            if dz2>r*r: continue
            mask=(xs[:,None]-x0)**2 + (ys[None,:]-y0)**2 + dz2 <= r*r
            plane|=mask
            if plane.all(): break
        occ+=int(plane.sum())
    return float(occ*vol)
def _vols(m,cid,probe=1.4,grid=0.5):
    conf=m.GetConformer(cid)
    coords=np.array([conf.GetAtomPosition(i) for i in range(m.GetNumAtoms())],float)
    radii=np.array([PT.GetRvdw(m.GetAtomWithIdx(i).GetAtomicNum()) for i in range(m.GetNumAtoms())],float)
    Vvdw=_grid_union(coords,radii,grid=grid); Vsas=_grid_union(coords,radii+probe,grid=grid); return Vvdw,Vsas
def _measure(name,smi,nconf=64,prune=0.5,T=298.15,grid=0.5,probe=1.4,seed=0x64C0FF):
    m=Chem.AddHs(Chem.MolFromSmiles(smi)); p=AllChem.ETKDGv3(); p.pruneRmsThresh=prune; p.randomSeed=seed
    ids=AllChem.EmbedMultipleConfs(m,numConfs=nconf,params=p); 
    if not ids: raise RuntimeError(f"Embed failed: {name}")
    opt=AllChem.MMFFOptimizeMoleculeConfs(m,maxIters=100); E=np.array([e for ok,e in opt],float)
    rows=[]
    for cid in ids:
        I1,I2,I3,NPR1,NPR2=_npr(m,cid); L,W,T=_extents(m,cid); Rg=_rg(m,cid,True); Vvdw,Vsas=_vols(m,cid,probe,grid)
        rows.append(dict(conf=int(cid),E=float(E[cid]),L=L,W=W,T=T,EndToEnd=_end2end(m,cid),Rg=Rg,NPR1=NPR1,NPR2=NPR2,V_vdw=Vvdw,V_sas_p1p4=Vsas))
    df=pd.DataFrame(rows).sort_values("E").reset_index(drop=True)
    kT=0.0019872041*T; w=np.exp(-(df.E-df.E.min())/kT); w=(w/w.sum()).to_numpy()
    def agg(col): a=df[col].to_numpy(); mu=float(np.average(a,weights=w)); sd=float(np.sqrt(np.average((a-mu)**2,weights=w))); return mu,sd
    out={"n_conformers":int(len(ids))}
    for k in ["L","W","T","EndToEnd","Rg","NPR1","NPR2","V_vdw","V_sas_p1p4"]:
        mu,sd=agg(k); out[f"{k}_mean"]=mu; out[f"{k}_sd"]=sd
    return out
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",default="FA SMILES.csv"); ap.add_argument("--out",default="FA_boltzmann_summaries.csv")
    ap.add_argument("--nconf",type=int,default=64); ap.add_argument("--grid",type=float,default=0.5)
    ap.add_argument("--probe",type=float,default=1.4); ap.add_argument("--temp",type=float,default=298.15)
    args=ap.parse_args()
    raw=pd.read_csv(args.input); cols={c.strip().lower():c for c in raw.columns}
    df_in=raw[[cols["fa"],cols["smiles"]]].rename(columns={cols["fa"]:"FA",cols["smiles"]:"SMILES"})
    df_in["FA"]=df_in["FA"].astype(str).str.strip(); df_in["SMILES"]=df_in["SMILES"].astype(str).str.strip()
    df_in=df_in[df_in["SMILES"]!=""].reset_index(drop=True)
    sums=[]
    for _,r in df_in.iterrows():
        s=_measure(r["FA"],r["SMILES"],nconf=args.nconf,grid=args.grid,probe=args.probe,T=args.temp); s["FA"]=r["FA"]; s["SMILES"]=r["SMILES"]; sums.append(s)
    cols_out=["FA","SMILES","n_conformers"]+[f"{k}_{stat}" for k in ["L","W","T","EndToEnd","Rg","NPR1","NPR2","V_vdw","V_sas_p1p4"] for stat in ["mean","sd"]]
    pd.DataFrame(sums,columns=cols_out).to_csv(args.out,index=False)
    meta=dict(rdkit_version=rdBase.rdkitVersion,python=sys.version.split()[0],platform=platform.platform(),
              date=str(datetime.date.today()),params=dict(n_conformers=args.nconf,grid_A=args.grid,probe_A=args.probe,temperature_K=args.temp))
    with open("run_metadata.json","w") as f: json.dump(meta,f,indent=2)
    print("Wrote:",args.out,"and run_metadata.json")
if __name__=="__main__": main()

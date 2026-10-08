import pandas as pd,numpy as np,json,re,collections
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier as HGC
from sklearn.model_selection import KFold,GroupKFold
from sklearn.metrics import roc_auc_score,brier_score_loss
u=pd.read_csv('/tmp/chk/up.tsv',sep='\t').fillna(''); u=u[u['Gene Names (primary)']!=''].drop_duplicates('Gene Names (primary)').reset_index(drop=True)
n=pd.read_csv('/tmp/het/nodes.tsv',sep='\t'); E=pd.read_csv('/tmp/het/edges.sif.gz',sep='\t'); gn=dict(zip(n.id,n.name)); tg=set(gn[x] for x in E[E.metaedge=='CbG'].target)
y=u['Gene Names (primary)'].isin(tg).astype(int).values
fam=u['Protein families'].str.split(',').str[0].str.strip().values; grp=np.array([f if f else 'single_'+a for f,a in zip(fam,u['Entry'])])
kw=u['Keywords'].str.split(';'); cnt=collections.Counter(k.strip() for ks in kw for k in ks if k.strip()); topk=[k for k,_ in cnt.most_common(300)]
KW=np.array([[k in set(x.strip() for x in ks) for k in topk] for ks in kw],float)
loc=u['Subcellular location [CC]'].str.replace('SUBCELLULAR LOCATION:','').str.replace(r'\{[^}]*\}','',regex=True).str.split(r'[.;,]'); lc=collections.Counter(t.strip().split('Note=')[0].strip() for ts in loc for t in ts if t.strip() and len(t.strip())<40); topl=[k for k,_ in lc.most_common(60)]
LC=np.array([[l in set(t.strip() for t in ts) for l in topl] for ts in loc],float)
aa='ACDEFGHIKLMNPQRSTVWY'; SQ=np.array([[s.count(a)/max(1,len(s)) for a in aa]+[np.log(len(s)+1)] for s in u['Sequence']])
FAMKW=set('Kinase,Receptor,G-protein coupled receptor,Ion channel,Transferase,Hydrolase,Protease,Oxidoreductase,Transducer,Voltage-gated channel,Ligand-gated ion channel,Serine/threonine-protein kinase,Tyrosine-protein kinase'.split(','))
KWab=KW[:,[i for i,k in enumerate(topk) if k not in FAMKW]]
R={'d1':dict(n=len(u),positives=int(y.sum()),rate=round(float(y.mean()),3),n_families=int(len(set(fam))-1),frac_with_family=round(float((fam!='').mean()),3),removed_keywords=[k for k in topk if k in FAMKW])}
def folds(kind,s):
    if kind=='random': return list(KFold(5,shuffle=True,random_state=s).split(KW))
    ug=sorted(set(grp)); r=np.random.default_rng(s); m=dict(zip(ug,r.permutation(len(ug)))); return list(GroupKFold(5).split(KW,groups=np.array([m[g] for g in grp])))
def run(X,kind,model='hgb',seeds=(0,1,2),shuffle=False,yy=None):
    yy=y if yy is None else yy; au=[]; oof=None
    for s in seeds:
        o=np.zeros(len(yy)); rng=np.random.default_rng(s)
        for tr,te in folds(kind,s):
            ytr=rng.permutation(yy[tr]) if shuffle else yy[tr]
            m=HGC(max_iter=150,max_depth=4,random_state=0) if model=='hgb' else LogisticRegression(C=.1,max_iter=500)
            m.fit(X[tr],ytr); o[te]=m.predict_proba(X[te])[:,1]
        au.append(roc_auc_score(yy,o)); oof=o
    return round(float(np.mean(au)),4),oof
XA=np.hstack([KW,LC,SQ]); XB=np.hstack([KWab,LC,SQ]); def_=lambda *a:None
def save(): json.dump(R,open('results/results.json','w'),indent=1,default=float)
R['d2_random']={'hgb':run(XA,'random')[0],'lr':run(XA,'random',model='lr')[0]}
# family prior baseline random
o=np.zeros(len(y)); 
for tr,te in folds('random',0):
    d=pd.DataFrame({'g':grp[tr],'y':y[tr]}).groupby('g').y.agg(['sum','count']); prior=y[tr].mean(); r=(d['sum']+prior*2)/(d['count']+2); o[te]=[r.get(g,prior) for g in grp[te]]
R['d2_random']['family_prior_baseline']=round(float(roc_auc_score(y,o)),4); print('d2',R['d2_random'],flush=True); save()
a,oo=run(XA,'cold_family'); al,_=run(XA,'cold_family',model='lr'); R['d3_cold_family']={'hgb':a,'lr':al,'best':max(a,al)}; R['d3_gate_ge_0.75']=bool(max(a,al)>=.75); print('d3',R['d3_cold_family'],flush=True); save()
a4,_=run(XB,'cold_family'); a4l,_=run(XB,'cold_family',model='lr'); R['d4_keyword_ablated_cold_family']={'hgb':a4,'lr':a4l,'best':max(a4,a4l)}; R['d4_gate_ge_0.70']=bool(max(a4,a4l)>=.70); print('d4',R['d4_keyword_ablated_cold_family'],flush=True); save()
R['d5_sequence_only_cold_family']={'lr':run(SQ,'cold_family',model='lr')[0],'hgb':run(SQ,'cold_family')[0]}
R['d6_groups_cold_family_hgb']={'keywords_only':run(KW,'cold_family')[0],'location_only':run(LC,'cold_family')[0],'sequence_only':R['d5_sequence_only_cold_family']['hgb']}; print('d5/6',R['d5_sequence_only_cold_family'],R['d6_groups_cold_family_hgb'],flush=True); save()
R['d7_shuffle_cold_family']=run(XA,'cold_family',seeds=(0,1,2),shuffle=True)[0]
cls={'kinase':'Kinase','receptor':'Receptor','ion_channel':'Ion channel','enzyme':'Hydrolase'}
R['d8_per_class_auroc_cold_family']={}
for k,w in cls.items():
    j=topk.index(w) if w in topk else None
    if j is None: continue
    mk=KW[:,j]==1; 
    if mk.sum()>50 and 0<y[mk].sum()<mk.sum(): R['d8_per_class_auroc_cold_family'][k]=dict(n=int(mk.sum()),pos=int(y[mk].sum()),auroc=round(float(roc_auc_score(y[mk],oo[mk])),3))
mk=~np.isin(np.arange(len(y)),np.where(KW[:,[topk.index(w) for w in cls.values() if w in topk]].sum(1)>0)[0]); R['d8_per_class_auroc_cold_family']['other']=dict(n=int(mk.sum()),pos=int(y[mk].sum()),auroc=round(float(roc_auc_score(y[mk],oo[mk])),3))
bs=brier_score_loss(y,oo); bins=np.linspace(0,1,11); ece=0
for i in range(10):
    m_=(oo>=bins[i])&(oo<bins[i+1]+(i==9)); 
    if m_.sum(): ece+=m_.sum()/len(y)*abs(oo[m_].mean()-y[m_].mean())
R['d9_calibration_cold_family']=dict(brier=round(float(bs),4),brier_baseline=round(float(y.mean()*(1-y.mean())),4),ece=round(float(ece),4))
un=np.where(y==0)[0]; R['d10_unlabeled_ranked']=dict(n_unlabeled=len(un),p_gt_0_5=int((oo[un]>.5).sum()),p_gt_0_8=int((oo[un]>.8).sum()),top25=[u['Gene Names (primary)'][i] for i in un[np.argsort(-oo[un])[:25]]],note='hypotheses from a label of old DrugBank binders; NOT validated druggable targets and no expansion claim')
save(); print(json.dumps(R,indent=1,default=float)[:3500])

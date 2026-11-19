import json, subprocess, sqlite3, os, pandas as pd
EXE=r"C:/Users/Administrator/AppData/Local/Programs/APSIM2025.3.7676.0/bin/Models.exe"
base=json.load(open('../apsim/SD_SF_A5B2.apsimx',encoding='utf-8-sig'))
def run(cmds, tag):
    j=json.loads(json.dumps(base))
    f=[x for x in j['Children'][0]['Children'] if x.get('Name')=='Field'][0]
    mz=[x for x in f['Children'] if x.get('Name')=='Maize'][0]
    mz['Children']=[{"$type":"Models.Core.Folder, Models","Name":"Cultivars","Children":[
        {"$type":"Models.PMF.Cultivar, Models","Command":list(cmds),"Name":"CN_TEST",
         "Children":[],"Enabled":True,"ReadOnly":False}],"Enabled":True,"ReadOnly":False}]
    for p in [c for c in f['Children'] if c.get('Name','').startswith('Sow')][0]['Parameters']:
        if p['Key']=='CultivarName': p['Value']='CN_TEST'
    fn=f'{tag}.apsimx'; json.dump(j,open(fn,'w',encoding='utf-8'),ensure_ascii=False)
    p=subprocess.run([EXE,fn],capture_output=True,text=True)
    err=p.stdout+p.stderr; db=fn.replace('.apsimx','.db')
    if 'ERROR' in err or not os.path.exists(db):
        return None, 'FAIL ' + str([l.strip() for l in err.splitlines() if 'Invalid path' in l or 'Unable' in l][:1])
    con=sqlite3.connect(db); df=pd.read_sql('select * from Report',con); con.close()
    st=df['Maize.Phenology.CurrentStageName'].astype(str)
    fl=df.loc[st.eq('Flowering'),'Clock.Today']
    mt=df.loc[st.eq('Maturity'),'Clock.Today']
    return dict(fl=str(fl.iloc[0])[:10] if len(fl) else None,
                mt=str(mt.iloc[0])[:10] if len(mt) else None,
                y=round(df['Maize.Grain.Total.Wt'].max()*10),
                lai=round(df['Maize.Leaf.LAI'].max(),2),
                bio=round(df['Maize.AboveGround.Wt'].max()*10)), 'OK'

cands = {
 "Phenology.Juvenile.Target.FixedValue": "=[340]",
 "Phenology.FlagLeafToFlowering.Target.FixedValue": "=[60]",
 "Phenology.FloweringToGrainFilling.Target.FixedValue": "=[140]",
 "Phenology.GrainFilling.Target.FixedValue": "=[600]",
 "Phenology.Maturing.Target.FixedValue": "=[15]",
 "Leaf.Photosynthesis.RUE.FixedValue": "=[2.4]",
 "Leaf.CohortParameters.MaxArea.MaximumArea.FixedValue": "=[0.0]",
 "Structure.FinalLeafNumber.ValueToHold.FixedValue": "=[20]",
 "Grain.MaximumPotentialGrainSize.FixedValue": "=[0.35]",
 "Grain.MaximumGrainsPerCob.FixedValue": "=[800]",
 "Leaf.CohortParameters.SenescenceAcceleration.One.FixedValue": "=[1.5]",
 "Leaf.CohortParameters.LagAcceleration.One.FixedValue": "=[1.2]",
 "Root.RootFrontVelocity.PotentialRootFrontVelocity.PreFlowering.RootFrontVelocity.FixedValue": "=[30]",
 "Leaf.FVPD.FixedValue": "=[0.95]",
}
ok=[]
for path, val in cands.items():
    expr = val.replace("=[", "=").replace("]", "")
    key = f"[Maize].{path} {val}"
    r, status = run([key], 'p'+str(abs(hash(path))%99999))
    print(f"{'OK ' if r else 'NO '} {path:72s} {status if not r else r}")
    if r: ok.append(path)
print("\nWRITABLE:", ok)

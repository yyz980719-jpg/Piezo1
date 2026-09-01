import sys, os, json, csv, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import net_helper as n

def safe_json(url, default):
    for _ in range(5):
        try:
            raw = n.ensembl_get(url)
            if raw.lstrip().startswith("{"):
                return json.loads(raw)
        except Exception:
            pass
        time.sleep(3)
    return default

def chembl(path):
    return safe_json("https://www.ebi.ac.uk/chembl/api/data/"+path,
                    {"page_meta":{"total_count":0},"targets":[],"mechanisms":[],"target_components":[]})

genes=["PIEZO1","BCL6","RELA","RXRB","NR1H2","GATAD2A","ATF6","TCF7L1","PKNOX2","ZNF853"]
rows=[]
for g in genes:
    try:
        time.sleep(1.0)
        t=chembl(f"target?target_synonym__iexact={g}")
        human=list(dict.fromkeys(x["target_chembl_id"] for x in t.get("targets",[]) if x.get("organism")=="Homo sapiens"))
        ttype=""; n_all=0; n4=0
        for cid in human:
            time.sleep(0.8)
            n_all+=chembl(f"mechanism?target_chembl_id={cid}&limit=1").get("page_meta",{}).get("total_count",0)
            time.sleep(0.8)
            n4+=chembl(f"mechanism?target_chembl_id={cid}&max_phase=4&limit=1").get("page_meta",{}).get("total_count",0)
            if not ttype:
                td=chembl(f"target/{cid}")
                ttype=td.get("target_type","")
        time.sleep(1.0)
        gwas=safe_json(f"https://www.ebi.ac.uk/gwas/rest/api/associations?geneName={g}&size=200",{"_embedded":{"associations":[]}})
        emb=gwas.get("_embedded",{}).get("associations",[])
        oa=0
        for a in emb:
            traits=[et.get("trait","") for et in a.get("efoTraits",[])]
            if a.get("traitName"): traits.append(a.get("traitName"))
            if any("osteoarthr" in (x or "").lower() for x in traits): oa+=1
        rows.append(dict(gene=g, n_human_targets=len(human), human_target_ids=";".join(human),
                         target_type=ttype, known_drug_mech=n_all, phase4_drugs=n4,
                         gwas_assoc_total=len(emb), gwas_oa_assoc=oa))
        print(f"{g:9s} nT={len(human)} type={ttype:20s} mech={n_all:5d} ph4={n4:4d} gwas={len(emb):4d} oa={oa}", flush=True)
    except Exception as e:
        print(f"{g} ERROR {e}", flush=True)
        rows.append(dict(gene=g, n_human_targets=-1, human_target_ids="", target_type="", known_drug_mech=-1, phase4_drugs=-1, gwas_assoc_total=-1, gwas_oa_assoc=-1))

os.makedirs("results", exist_ok=True)
with open("results/29_s16_chembl_verify.csv","w",newline="",encoding="utf-8-sig") as fh:
    w=csv.DictWriter(fh, fieldnames=["gene","n_human_targets","human_target_ids","target_type","known_drug_mech","phase4_drugs","gwas_assoc_total","gwas_oa_assoc"])
    w.writeheader(); w.writerows(rows)
print("WROTE results/29_s16_chembl_verify.csv", flush=True)

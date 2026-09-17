"""Bounded V63 archive workflow. Never prints credentials or private share links."""
import os,json,time,hashlib,zipfile,io
from datetime import datetime,timezone
import requests
ROOT='https://zenodo.org/api'
OLD='22231510'
NAME='PIEZO1_V63_public_reproducibility.zip'
SHA='be39084d9bdfaf7e46244d71fb873d862fbb43d7017927be02e4e4347c8bba82'
RELEASE='https://github.com/yyz980719-jpg/Piezo1/releases/tag/v63-reproducibility'
S=requests.Session()
S.headers['Authorization']='Bearer '+os.environ['ZENODO_TOKEN']
def req(method,url,**kw):
    assert url.startswith(ROOT+'/')
    for attempt in range(3 if method=='GET' else 1):
        r=S.request(method,url,timeout=(20,180),allow_redirects=False,**kw)
        if r.status_code not in {429,502,503,504}:break
        if attempt<2:time.sleep(5*(attempt+1))
    if not r.ok:raise RuntimeError(f'Zenodo {method} returned HTTP {r.status_code}')
    return r.json() if r.content else None
def brief(d):
    m=d.get('metadata',{})
    return dict(id=d.get('id'),state=d.get('state'),submitted=d.get('submitted'),doi=d.get('doi',m.get('doi')),conceptrecid=d.get('conceptrecid'),title=m.get('title'),version=m.get('version'),creators=m.get('creators'),license=m.get('license'),files=[dict(name=f.get('filename',f.get('key')),id=f.get('id'),bytes=f.get('filesize',f.get('size')),checksum=f.get('checksum')) for f in d.get('files',[])])
def asset():
    r=requests.get(RELEASE.replace('/tag/','/download/')+'/'+NAME,timeout=(20,180))
    r.raise_for_status()
    data=r.content
    assert len(data)==79417872 and hashlib.sha256(data).hexdigest()==SHA
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        assert z.testzip() is None
        assert not any(n.lower().endswith(('.docx','.doc')) for n in z.namelist())
    print('Release asset verified:',len(data),'bytes; sha256',SHA,flush=True)
    return data

def metadata(old):
    return dict(title='PIEZO1 and TRPV4 transcript associations in osteoarthritic cartilage: V63 selected-cell reproducibility materials',
      upload_type='software',version='v63-reproducibility',access_right='open',license=old['metadata']['license'],
      creators=old['metadata']['creators'],publication_date=datetime.now(timezone.utc).date().isoformat(),
      description='<p>Version-fixed code, selected-cell input arrays, donor and sample mappings, numerical outputs, figure source data, direct dependency versions and replay instructions for the V63 PIEZO1-TRPV4 cartilage association analyses. The subsequent V64 manuscript revision is editorial and uses these same analyses.</p><p>The inputs derive from GEO GSE255460 and GSE220243. Replays cover selected-cell correlations, donor summaries, fixed statistical families, measurement sensitivity, member and inherited-state comparisons, and nested cell-resampling diagnostics. They do not re-run raw-read alignment, initial cell calling, upstream annotation or doublet-model fitting. Internal computational checks are not independent biological validation.</p><p>The main and supplementary manuscripts are not included. The MIT license applies to original code; third-party source data and annotations retain their original rights and attribution, as detailed in THIRD_PARTY_NOTICES.md.</p><p>The named ZIP is byte-identical to the published GitHub V63 release asset. Earlier versions remain separate historical records.</p>',
      keywords=['PIEZO1','TRPV4','osteoarthritis','cartilage','single-cell RNA sequencing','donor-level inference','reproducible research'],
      related_identifiers=[dict(identifier=RELEASE,relation='isSupplementTo',scheme='url')]+[dict(identifier='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+a,relation='isDerivedFrom',scheme='url') for a in ['GSE255460','GSE220243']],
      notes='Code creator credited here; manuscript authorship is separate. Fixed analysis commit: 14ab064c5a58df5fe6b07a223a07fa6897ce136e. Archive SHA-256: '+SHA+'.')

def validate(d,data):
    assert str(d['id'])!=OLD and str(d['conceptrecid'])=='22231509'
    assert d['metadata']['version']=='v63-reproducibility'
    assert d['metadata']['access_right']=='open'
    assert 'manuscripts are not included' in d['metadata']['description']
    assert len(d['files'])==1
    f=d['files'][0]
    assert f.get('filename',f.get('key'))==NAME
    assert int(f.get('filesize',f.get('size')))==len(data)
    assert f['checksum'].removeprefix('md5:')==hashlib.md5(data).hexdigest()

def main():
    action=os.environ.get('ARCHIVE_ACTION','inspect')
    d=req('GET',f'{ROOT}/deposit/depositions/{OLD}')
    print(json.dumps(brief(d),indent=2))
    for key in ['latest','latest_draft']:
        link=d.get('links',{}).get(key)
        if link:
            r=req('GET',link)
            print(key,json.dumps(brief(r),indent=2))
    if action=='inspect':return
    assert action in {'prepare','publish'}
    data=asset()
    if action=='prepare':
        latest=req('GET',d['links']['latest'])
        assert str(latest['id'])==OLD, 'Newer version exists; inspect before mutation'
        draft=req('GET',d['links']['latest_draft'])
        if str(draft['id'])==OLD:
            result=req('POST',f'{ROOT}/deposit/depositions/{OLD}/actions/newversion')
            draft=req('GET',result['links']['latest_draft'])
            print('Created new draft',draft['id'],flush=True)
        else:
            assert draft['metadata'].get('version')=='v63-reproducibility', 'Unknown draft: manual review required'
        assert not draft['submitted'] and str(draft['id'])!=OLD
        assert str(draft['conceptrecid'])=='22231509'
        target=f"{ROOT}/deposit/depositions/{draft['id']}"
        draft=req('PUT',target,json={'metadata':metadata(d)})
        # Delete only inherited snapshots in this explicitly verified unpublished new version.
        inherited={f['id'] for f in d['files']}
        for f in draft.get('files',[]):
            name=f.get('filename',f.get('key'))
            if name==NAME:continue
            assert f['id'] in inherited, 'Unexpected draft file; do not delete'
            req('DELETE',target+'/files/'+f['id'])
            print('Removed inherited draft snapshot:',name,flush=True)
        draft=req('GET',target)
        if not any(f.get('filename',f.get('key'))==NAME for f in draft['files']):
            req('PUT',draft['links']['bucket']+'/'+NAME,data=data,headers={'Content-Type':'application/octet-stream'})
        draft=req('GET',target)
        validate(draft,data)
        print('VERIFIED_DRAFT',json.dumps(brief(draft),indent=2),flush=True)
    else:
        draft_id=os.environ['DRAFT_ID']
        assert draft_id.isdigit() and draft_id!=OLD
        target=f'{ROOT}/deposit/depositions/{draft_id}'
        draft=req('GET',target)
        validate(draft,data)
        assert draft['metadata']['creators']==d['metadata']['creators']
        if not draft['submitted']:
            draft=req('POST',target+'/actions/publish')
        validate(draft,data)
        assert draft['submitted']
        print('PUBLISHED',json.dumps(brief(draft),indent=2),flush=True)
if __name__=='__main__':main()

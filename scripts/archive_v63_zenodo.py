"""Bounded V63 archive workflow. Never prints credentials or private share links."""
import os,json,time
import requests
ROOT='https://zenodo.org/api'
OLD='22231510'
S=requests.Session()
S.headers['Authorization']='Bearer '+os.environ['ZENODO_TOKEN']
def req(method,url,**kw):
    assert url.startswith(ROOT+'/')
    for attempt in range(3 if method=='GET' else 1):
        r=S.request(method,url,timeout=(20,60),**kw)
        if r.status_code not in {429,502,503,504}:break
        if attempt<2:time.sleep(5*(attempt+1))
    if not r.ok:raise RuntimeError(f'Zenodo {method} returned HTTP {r.status_code}')
    return r.json() if r.content else None
def brief(d):
    m=d.get('metadata',{})
    return dict(id=d.get('id'),state=d.get('state'),submitted=d.get('submitted'),doi=d.get('doi',m.get('doi')),conceptrecid=d.get('conceptrecid'),title=m.get('title'),version=m.get('version'),creators=m.get('creators'),license=m.get('license'),files=[dict(name=f.get('filename',f.get('key')),id=f.get('id'),bytes=f.get('filesize',f.get('size')),checksum=f.get('checksum')) for f in d.get('files',[])])
def main():
    d=req('GET',f'{ROOT}/deposit/depositions/{OLD}')
    print(json.dumps(brief(d),indent=2))
    for key in ['latest','latest_draft']:
        link=d.get('links',{}).get(key)
        if link:
            r=req('GET',link)
            print(key,json.dumps(brief(r),indent=2))
    assert os.environ.get('ARCHIVE_ACTION','inspect')=='inspect'
if __name__=='__main__':main()

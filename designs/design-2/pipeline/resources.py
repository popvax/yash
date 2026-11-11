"""Prepare verified public resources without changing the original reference."""
from pathlib import Path
import base64
import hashlib
import json
import sys
import urllib.request
import urllib.parse

import pysam
from common import sha,save,receipt,now

ROOT=Path(sys.argv[1]);OUT=ROOT/"resources";OUT.mkdir(parents=True,exist_ok=True)
BASE="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/"
NAMES=["Homo_sapiens_assembly38.haplotype_database.txt",
       "Homo_sapiens_assembly38.dbsnp138.vcf.gz",
       "Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi",
       "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
       "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"]


def main():
    records=[]
    for name in NAMES:
        target=OUT/name;done=OUT/(name+".receipt.json")
        if done.exists():
            record=json.loads(done.read_text())
            assert target.stat().st_size==record["bytes"] and sha(target)==record["sha256"]
            records.append(record);continue
        url=BASE+name
        metadata_url=("https://storage.googleapis.com/storage/v1/b/gcp-public-data--broad-references/o/"
                      +urllib.parse.quote("hg38/v0/"+name,safe="")+"?fields=md5Hash,size,generation,contentEncoding")
        with urllib.request.urlopen(metadata_url,timeout=60) as response:
            metadata=json.load(response)
        pinned_url=url+"?generation="+metadata["generation"]
        save(OUT/"status.json",dict(phase="downloading",file=name,updated_utc=now()))
        partial=target.with_name(target.name+".part")
        digest=hashlib.sha256();md5=hashlib.md5();total=0
        encoding="gzip" if metadata.get("contentEncoding")=="gzip" else "identity"
        request=urllib.request.Request(pinned_url,headers={"Accept-Encoding":encoding})
        with urllib.request.urlopen(request,timeout=120) as response, partial.open("wb") as stream:
            headers=dict(response.headers.items());expected=int(metadata["size"])
            for block in iter(lambda:response.read(8*1024*1024),b""):
                stream.write(block);digest.update(block);md5.update(block);total+=len(block)
        assert total==expected
        actual_md5=base64.b64encode(md5.digest()).decode()
        assert metadata.get("md5Hash")==actual_md5,"Published resource MD5 is unavailable or mismatched"
        partial.rename(target)
        record=dict(name=name,url=url,bytes=total,sha256=digest.hexdigest(),
                    md5_base64=actual_md5,provider_md5_verified=True,http_headers=headers,
                    provider_metadata=metadata,metadata_url=metadata_url,pinned_download_url=pinned_url)
        save(done,record);records.append(record)
        print(json.dumps({"resource":name,"bytes":total,"verified":True}),flush=True)
    # These VCFs are BQSR masks, not somatic-calling nomination resources.
    for name in [NAMES[1],NAMES[3]]:
        with pysam.VariantFile(OUT/name) as variants:
            next(variants.fetch("chr1",10000,100000),None)
    receipt(OUT,[OUT/name for name in NAMES],resources=records,
            scope="Broad hg38 haplotype fingerprint map and known-variant BQSR masks",
            source_reference_unchanged=True,worker_sha256=sha(__file__))


if __name__=="__main__":main()

"""Publish an already verified Design #2 transfer into the existing data catalog.

Requires full-GET transfer verification. Preserves prior objects and bucket
permissions; changes INDEX.json conditionally after immutable metadata is checked.
"""
from pathlib import Path
import argparse
import concurrent.futures
import datetime
import gzip
import hashlib
import json

from common import load_credentials, client, get_bytes
from botocore.exceptions import ClientError

REPO = Path(__file__).resolve().parents[2]
PREFIX = "metadata/2026-11-07/design-2/"


def encoded(value):
    return (json.dumps(value, indent=2) + "\n").encode()


def sha(body):
    return hashlib.sha256(body).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, default=REPO / "work/design-2-r2")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    plan = json.loads((args.work / "plan.json").read_text())
    verified = json.loads((args.work / "verification.json").read_text())
    assert verified["status"] == "VERIFIED_BY_FULL_GET_SHA256"
    assert verified["plan_sha256"] == sha((args.work / "plan.json").read_bytes())
    assert not verified["original_normal_FASTQs_included"]
    assert set(verified["objects"]) == {o["key"] for o in plan["objects"]}
    for item in plan["objects"]:
        proof = verified["objects"][item["key"]]
        assert proof["sha256"] == item["sha256"] and proof["bytes"] == item["bytes"]
        assert proof["verification"] == "FULL_GET_SHA256" and proof["storage_class"] == "STANDARD"
    baseline = json.loads((args.work / "before-bundled-manifest.json").read_text())
    before = json.loads((args.work / "before-index.json").read_text())
    previous_remote = json.loads((args.work / "before-remote-manifest.json").read_text())
    previous_objects = {r["key"]: r for r in previous_remote["objects"]}
    assert all(previous_objects[r["key"]]["sha256"] == r["sha256"] and previous_objects[r["key"]]["bytes"] == r["bytes"] for r in baseline["objects"])
    objects = baseline["objects"] + plan["objects"]
    assert len({r["key"] for r in objects}) == len(objects)
    manifest = dict(baseline, snapshot_date="2026-11-07", objects=objects,
                    scope="Recorded tumour/RNA inputs and references, plus matched-normal alignment and Design #2 outputs; original normal FASTQs excluded")
    local_current = json.loads((REPO / "data/catalog/object-manifest.json").read_text())
    assert local_current in [baseline, manifest], "Local catalog changed concurrently; preserve the other edit"
    blob = encoded(manifest)
    print(json.dumps(dict(execute=args.execute, catalog_objects=len(objects), new_objects=len(plan["objects"]),
                          new_bytes=sum(r["bytes"] for r in plan["objects"]), bucket_permissions_changed=False), indent=2))
    if not args.execute:
        return
    load_credentials(); s3 = client(); bucket = plan["bucket"]
    current = s3.get_object(Bucket=bucket, Key="INDEX.json")
    try:
        current_index = json.loads(current["Body"].read())
    finally:
        current["Body"].close()
    assert current_index == before or current_index.get("object_manifest_sha256") == sha(blob), "Remote index changed concurrently"

    def check(item):
        metadata = s3.head_object(Bucket=bucket, Key=item["key"])
        assert metadata["ContentLength"] == item["bytes"] and metadata["ETag"] == verified["objects"][item["key"]]["ETag"], item["key"]
        assert metadata.get("StorageClass", "STANDARD") == "STANDARD"
    with concurrent.futures.ThreadPoolExecutor(max_workers=12) as pool:
        list(pool.map(check, plan["objects"]))

    source = b"".join((json.dumps(dict(path=o["source_path"], type="1", bytes=o["bytes"], storage="individual", key=o["key"], sha256=o["sha256"], uri=o["uri"])) + "\n").encode() for o in plan["objects"])
    source = gzip.compress(source, mtime=0)
    publication = dict(status="VERIFIED_BY_FULL_GET_SHA256", new_objects=verified["objects"],
        verified_manifest_objects=len(objects), new_object_count=len(plan["objects"]), inherited_object_count=len(baseline["objects"]),
        inherited_verification_receipt_key=before["verification_receipt_key"],
        inherited_verification_receipt_sha256=before["verification_receipt_sha256"],
        inherited_remote_manifest_key=before["object_manifest_key"], inherited_remote_manifest_sha256=before["object_manifest_sha256"],
        source_transfer_receipts=verified.get("source_transfer_receipts", []),
        bundled_manifest_sha256=sha(blob), bucket_permissions_changed=False, original_normal_FASTQs_included=False,
        historical_objects_deleted=False, verified_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    prepared = args.work / "prepared-publication.json"
    if prepared.exists():
        recorded = json.loads(prepared.read_text())
        assert recorded["new_objects"] == publication["new_objects"] and recorded["bundled_manifest_sha256"] == sha(blob)
        publication = recorded
    else:
        prepared.write_bytes(encoded(publication))
    metadata = {"object-manifest.json": blob, "source-manifest.jsonl.gz": source,
                "verification.json": encoded(publication), "previous-INDEX.json": encoded(before)}
    for name, body in metadata.items():
        key = PREFIX + name
        try:
            s3.put_object(Bucket=bucket, Key=key, Body=body, ContentType="application/gzip" if name.endswith(".gz") else "application/json", IfNoneMatch="*")
        except ClientError as error:
            assert error.response["Error"]["Code"] in ["PreconditionFailed", "412"], error
        assert get_bytes(s3, key) == body, "Published metadata mismatch: " + key

    index = dict(before, snapshot_date="2026-11-07", lead_design="design-2",
        object_manifest_key=PREFIX + "object-manifest.json", object_manifest_sha256=sha(blob),
        verification_receipt_key=PREFIX + "verification.json", verification_receipt_sha256=sha(metadata["verification.json"]),
        verification_status="VERIFIED_BY_FULL_GET_SHA256", verified_manifest_objects=len(objects), verified_utc=publication["verified_utc"],
        additional_source_manifests=[dict(key=PREFIX + "source-manifest.jsonl.gz", sha256=sha(source), scope="Matched-normal and Design #2 outputs")],
        previous_index_key=PREFIX + "previous-INDEX.json", original_normal_FASTQs_included=False,
        repository_snapshot_note="The 2026-10-31 repository archive is historical; current design documentation is in the GitHub repository")
    assert current_index == before or current_index == index, "Index fields changed concurrently; preserve the other publisher's update"
    s3.put_object(Bucket=bucket, Key="INDEX.json", Body=encoded(index), ContentType="application/json",
                  CacheControl="no-store", IfMatch=current["ETag"])
    assert json.loads(get_bytes(s3, "INDEX.json")) == index
    assert sha(get_bytes(s3, index["object_manifest_key"])) == index["object_manifest_sha256"]

    (REPO / "data/catalog/object-manifest.json").write_bytes(blob)
    access_path = REPO / "data/catalog/access.json"; access = json.loads(access_path.read_text())
    if access.get("bundled_manifest", "object-manifest.json") == "object-manifest.json":
        access["bundled_manifest_sha256"] = sha(blob)
    else:
        # Public hosting is now curated/deployed separately from private archives.
        access["native_manifest"] = "object-manifest.json"
        access["native_manifest_sha256"] = sha(blob)
    access_path.write_bytes(encoded(access))
    source_path = REPO / "data/catalog/sources.json"; sources = json.loads(source_path.read_text())
    added = [dict(key=o["key"], uri=o["uri"], sha256=o["sha256"], bytes=o["bytes"]) for o in plan["objects"] if o["category"] == "raw"]
    sources["inputs"] = [r for r in sources["inputs"] if r["key"] not in {o["key"] for o in added}] + added
    sources["matched_normal_and_design2_manifest"] = "../matched-normal/manifest.json"
    source_path.write_bytes(encoded(sources))
    target = REPO / "data/matched-normal"; target.mkdir(exist_ok=True)
    summary = dict(schema_version=1, case="TN25-112103", normal_sample="TN25-112103-N", design="design-2",
                   source_type="White-blood-cell genomic DNA", sequencing_instrument="Unknown", reference_sha256="c3e61ad6d770db22a6f27fdbcd3dba93d0d80f649018c8887ec803f12ea90450",
                   normal_genotyping_scope="Primary nuclear annotated CDS/start-stop intervals plus 100-bp flanks; not whole-genome genotyping",
                   original_normal_FASTQs_included=False, bucket_access_policy_changed=False, objects=plan["objects"])
    (target / "manifest.json").write_bytes(encoded(summary))
    (target / "verification.json").write_bytes(encoded(publication))
    (args.work / "published-index.json").write_bytes(encoded(index))
    print("PUBLISHED_AND_VERIFIED", len(objects), "catalog objects; bucket permissions unchanged")


if __name__ == "__main__":
    main()

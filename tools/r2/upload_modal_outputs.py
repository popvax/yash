"""Upload an explicit output plan from a Modal volume; verify full R2 GET hashes.

Credentials are provided only through the sandbox environment. Multipart state,
progress and receipts contain no credentials. Existing different objects are
never overwritten. This worker does not publish a catalog or change permissions.
"""
from pathlib import Path
import base64
import concurrent.futures
import hashlib
import json
import math
import os
import sys
import time

import boto3
from botocore.config import Config
from botocore.exceptions import ClientError
from ranged_readback import RangedBody

PART = 64 * 1024 * 1024
WORKERS = 8


def save(path, value):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2) + "\n")
    tmp.replace(path)


def main():
    plan = json.loads(Path(sys.argv[1]).read_text())
    out = Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
    assert plan["bucket"] == os.environ["R2_BUCKET"] == "yash-rna-foundry"
    assert plan["original_normal_FASTQs_included"] is False
    client = boto3.client("s3", endpoint_url=os.environ["R2_ENDPOINT"], region_name="auto",
        aws_access_key_id=os.environ["R2_ACCESS_KEY_ID"], aws_secret_access_key=os.environ["R2_SECRET_ACCESS_KEY"],
        config=Config(signature_version="s3v4", max_pool_connections=WORKERS + 4,
            retries={"mode": "adaptive", "max_attempts": 10},
            request_checksum_calculation="when_required", response_checksum_validation="when_required"))
    bucket = plan["bucket"]
    started = time.time(); total = sum(r["bytes"] for r in plan["objects"])
    receipt_path = out / "verification.json"
    receipt = json.loads(receipt_path.read_text()) if receipt_path.exists() else dict(objects={})
    last = 0

    def progress(item, phase, amount, force=False):
        nonlocal last
        if force or time.time() - last >= 5:
            last = time.time()
            value = dict(phase=phase, key=item["key"], object_bytes=item["bytes"], processed_bytes=amount,
                         objects_verified=len(receipt["objects"]), total_objects=len(plan["objects"]),
                         total_bytes=total, elapsed_seconds=round(time.time() - started, 1), updated_unix=time.time())
            save(out / "progress.json", value)
            if force or int(value["elapsed_seconds"]) % 30 < 5:
                print(json.dumps(value), flush=True)

    def head(key):
        try:
            return client.head_object(Bucket=bucket, Key=key)
        except ClientError as error:
            if error.response["Error"]["Code"] in ["404", "NoSuchKey", "NotFound"]:
                return None
            raise

    def verify(item, metadata):
        assert metadata["ContentLength"] == item["bytes"], item["key"]
        prior = receipt["objects"].get(item["key"])
        if prior and prior["bytes"] == item["bytes"] and prior["sha256"] == item["sha256"] and prior["ETag"] == metadata["ETag"]:
            return prior
        digest = hashlib.sha256(); size = 0
        progress(item, "VERIFY_FULL_R2_GET", 0, True)
        ranged = item["bytes"] >= 256*1024*1024
        body = (RangedBody(client, bucket, item["key"], item["bytes"], metadata["ETag"], workers=WORKERS)
                if ranged else client.get_object(Bucket=bucket, Key=item["key"])["Body"])
        try:
            for block in body.iter_chunks(chunk_size=8 * 1024 * 1024):
                digest.update(block); size += len(block); progress(item, "VERIFY_FULL_R2_GET", size)
        finally:
            body.close()
        assert size == item["bytes"] and digest.hexdigest() == item["sha256"], "R2 full-GET checksum mismatch: " + item["key"]
        return dict(bytes=size, sha256=digest.hexdigest(), ETag=metadata["ETag"],
                    storage_class=metadata.get("StorageClass", "STANDARD"), verification="FULL_GET_SHA256",
                    readback_mode="ORDERED_ETAG_PINNED_RANGE_GETS" if ranged else "SINGLE_GET", verified_unix=time.time())

    for item in plan["objects"]:
        assert item["key"].startswith(("raw/TN25-112103/normal/", "analysis/TN25-112103/design-2/"))
        assert not any(item["key"].endswith(s) for s in [".fastq.gz", ".fq.gz", ".fastq", ".fq"])
        path = Path(item["source_path"])
        assert str(path).startswith("/data/results/TN25-112103/matched-normal-20261107/")
        assert path.stat().st_size == item["bytes"], item["source_path"]
        existing = head(item["key"])
        if existing is None and item["bytes"] < PART:
            data = path.read_bytes()
            assert hashlib.sha256(data).hexdigest() == item["sha256"], "Source checksum mismatch"
            progress(item, "UPLOAD_AND_SOURCE_SHA256", len(data), True)
            client.put_object(Bucket=bucket, Key=item["key"], Body=data, StorageClass="STANDARD",
                Metadata={"sha256": item["sha256"]}, ContentType="application/octet-stream", IfNoneMatch="*",
                ContentMD5=base64.b64encode(hashlib.md5(data).digest()).decode("ascii"))
            existing = head(item["key"])
        if existing is None:
            state_path = out / (hashlib.sha256(item["key"].encode()).hexdigest() + ".multipart.json")
            state = json.loads(state_path.read_text()) if state_path.exists() else None
            if state:
                assert state["key"] == item["key"] and state["sha256"] == item["sha256"]
            else:
                upload = client.create_multipart_upload(Bucket=bucket, Key=item["key"], StorageClass="STANDARD",
                    Metadata={"sha256": item["sha256"]}, ContentType="application/octet-stream")
                state = dict(key=item["key"], sha256=item["sha256"], upload_id=upload["UploadId"])
                save(state_path, state)
            remote_parts = {}
            paginator = client.get_paginator("list_parts")
            for page in paginator.paginate(Bucket=bucket, Key=item["key"], UploadId=state["upload_id"]):
                for part in page.get("Parts", []):
                    remote_parts[part["PartNumber"]] = part
            completed = {}; pending = {}; digest = hashlib.sha256(); processed = 0

            def put(number, block):
                md5 = base64.b64encode(hashlib.md5(block).digest()).decode("ascii")
                result = client.upload_part(Bucket=bucket, Key=item["key"], UploadId=state["upload_id"],
                    PartNumber=number, Body=block, ContentMD5=md5)
                return dict(PartNumber=number, ETag=result["ETag"])

            def collect(done):
                for future in done:
                    completed[pending.pop(future)] = future.result()

            progress(item, "UPLOAD_AND_SOURCE_SHA256", 0, True)
            with path.open("rb") as source, concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as pool:
                number = 0
                while block := source.read(PART):
                    number += 1; digest.update(block); processed += len(block)
                    prior = remote_parts.get(number)
                    if prior:
                        assert prior["Size"] == len(block)
                        assert prior["ETag"].strip('"') == hashlib.md5(block).hexdigest(), "Resume part differs from source"
                        completed[number] = dict(PartNumber=number, ETag=prior["ETag"])
                    else:
                        pending[pool.submit(put, number, block)] = number
                    if len(pending) >= WORKERS:
                        done, _ = concurrent.futures.wait(pending, return_when=concurrent.futures.FIRST_COMPLETED)
                        collect(done)
                    progress(item, "UPLOAD_AND_SOURCE_SHA256", processed)
                collect(list(pending))
            assert processed == item["bytes"] and digest.hexdigest() == item["sha256"], "Source checksum mismatch"
            assert len(completed) == math.ceil(item["bytes"] / PART)
            client.complete_multipart_upload(Bucket=bucket, Key=item["key"], UploadId=state["upload_id"],
                MultipartUpload={"Parts": [completed[n] for n in sorted(completed)]}, IfNoneMatch="*")
            state_path.unlink()
            existing = head(item["key"])
        receipt["objects"][item["key"]] = verify(item, existing)
        save(receipt_path, receipt)
        progress(item, "OBJECT_VERIFIED", item["bytes"], True)
    receipt.update(status="VERIFIED_BY_FULL_GET_SHA256", total_objects=len(plan["objects"]), total_bytes=total,
                   original_normal_FASTQs_included=False, storage_class="STANDARD", bucket_permissions_changed=False,
                   source_files_deleted=False, plan_sha256=hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest())
    save(receipt_path, receipt)
    print(json.dumps({k:v for k,v in receipt.items() if k != "objects"}), flush=True)


if __name__ == "__main__":
    main()

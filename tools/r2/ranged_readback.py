"""Bounded parallel S3 range reads, yielded in byte order for a full-object hash."""
from concurrent.futures import ThreadPoolExecutor


class RangedBody:
    def __init__(self, client, bucket, key, size, etag, workers=8, part_bytes=64*1024*1024):
        self.client, self.bucket, self.key = client, bucket, key
        self.size, self.etag = size, etag
        self.workers, self.part_bytes = workers, part_bytes
        self.pool = None

    def fetch(self, number):
        start = number * self.part_bytes; end = min(self.size, start+self.part_bytes)-1
        response = self.client.get_object(Bucket=self.bucket, Key=self.key,
                                         Range=f"bytes={start}-{end}", IfMatch=self.etag)
        body = response["Body"]
        try:
            assert response["ContentRange"] == f"bytes {start}-{end}/{self.size}"
            assert response["ETag"] == self.etag
            block = body.read()
        finally:
            body.close()
        assert len(block) == end-start+1
        return block

    def iter_chunks(self, chunk_size):
        count = (self.size+self.part_bytes-1)//self.part_bytes
        self.pool = ThreadPoolExecutor(max_workers=self.workers)
        pending = {}; next_number = 0
        for number in range(min(self.workers, count)):
            pending[number] = self.pool.submit(self.fetch, number); next_number += 1
        for number in range(count):
            block = pending.pop(number).result()
            # Refill only after the previous ordered block is consumed. Memory
            # remains bounded; requests pin the same object version with ETag.
            yield block
            if next_number < count:
                pending[next_number] = self.pool.submit(self.fetch, next_number); next_number += 1

    def close(self):
        if self.pool:
            self.pool.shutdown(wait=True, cancel_futures=True)
            self.pool = None

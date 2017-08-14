import hashlib
import json


def default_hash_fn(bytes_):
    return hashlib.blake2b(bytes_, digest_size=16)


def hash_obj(obj=None, hash_fn=default_hash_fn):
    return hash_fn(bytify_obj(obj=obj)).hexdigest()


def bytify_obj(obj=None):
    return json.dumps(obj, sort_keys=True).encode()

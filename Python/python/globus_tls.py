# SAGEhydrology | Python source
# SITEhydrology | Python source
# File: python/globus_tls.py
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
(
    'Use Windows-trusted public CA certificates for Globus HTTPS '
    'verification.'
)

import os, sys, ssl, hashlib, tempfile
from pathlib import Path


def environment():
    env = dict(os.environ)
    if sys.platform != "win32" or env.get("GLOBUS_SDK_VERIFY_SSL"):
        return env
    if env.get("REQUESTS_CA_BUNDLE"):
        env["GLOBUS_SDK_VERIFY_SSL"] = env["REQUESTS_CA_BUNDLE"]
        return env
    import certifi

    pieces = [Path(certifi.where()).read_text(encoding="ascii")]
    for store in ("ROOT", "CA"):
        for certificate, encoding, trust in ssl.enum_certificates(store):
            if encoding == "x509_asn" and (
                trust is True or "1.3.6.1.5.5.7.3.1" in trust
            ):
                pieces.append(ssl.DER_cert_to_PEM_cert(certificate))
    content = "\n".join(pieces)
    folder = (
        Path(os.environ.get("LOCALAPPDATA", tempfile.gettempdir()))
        / "SAGEhydrology/Python/trusted_ca"
    )
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / (
        hashlib.sha256(content.encode("ascii")).hexdigest() + ".pem"
    )
    if not path.exists():
        path.write_text(content, encoding="ascii")
    env["GLOBUS_SDK_VERIFY_SSL"] = str(path)
    env["REQUESTS_CA_BUNDLE"] = str(path)
    return env

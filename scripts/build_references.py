#!/usr/bin/env python3

from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import datetime as dt
import json
import os
import re
import ssl
import sys
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any


DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")
BIB_ENTRY_RE = re.compile(r"^@(\w+)\s*{\s*([^,]+)\s*,", flags=re.IGNORECASE)


def utc_now_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(r) for r in reader]


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


@dataclass(frozen=True)
class FetchResult:
    ok: bool
    status: int | None
    url: str
    final_url: str
    content_type: str
    body: bytes
    error: str


def fetch_url(
    url: str,
    headers: dict[str, str] | None = None,
    timeout_sec: float = 20.0,
) -> FetchResult:
    ssl_context = None
    try:
        cafile_candidates = [
            os.environ.get("SSL_CERT_FILE", ""),
            "/etc/ssl/cert.pem",  # macOS (commonly present)
            "/etc/ssl/certs/ca-certificates.crt",  # Debian/Ubuntu
            "/etc/pki/tls/certs/ca-bundle.crt",  # RHEL/CentOS/Fedora
        ]
        for cafile in cafile_candidates:
            if cafile and Path(cafile).exists():
                ssl_context = ssl.create_default_context(cafile=cafile)
                break
    except Exception:
        ssl_context = None

    req = urllib.request.Request(url, headers=headers or {})
    try:
        with urllib.request.urlopen(req, timeout=timeout_sec, context=ssl_context) as resp:
            body = resp.read()
            status = getattr(resp, "status", None)
            final_url = getattr(resp, "url", url)
            content_type = (resp.headers.get("Content-Type") or "").strip()
            return FetchResult(
                ok=True,
                status=status,
                url=url,
                final_url=final_url,
                content_type=content_type,
                body=body,
                error="",
            )
    except urllib.error.HTTPError as e:
        return FetchResult(
            ok=False,
            status=e.code,
            url=url,
            final_url=url,
            content_type=(e.headers.get("Content-Type") or "").strip() if e.headers else "",
            body=e.read() if hasattr(e, "read") else b"",
            error=f"HTTPError {e.code}",
        )
    except Exception as e:  # noqa: BLE001 - tooling script
        return FetchResult(
            ok=False,
            status=None,
            url=url,
            final_url=url,
            content_type="",
            body=b"",
            error=f"{type(e).__name__}: {e}",
        )


def normalize_doi(doi: str) -> str:
    doi = (doi or "").strip()
    doi = doi.replace("https://doi.org/", "").replace("http://doi.org/", "")
    return doi


def rewrite_bibtex_key(bibtex: str, citekey: str) -> str:
    bibtex = bibtex.lstrip("\ufeff").strip()
    match = BIB_ENTRY_RE.search(bibtex)
    if not match:
        return bibtex
    entry_type = match.group(1)
    # Replace only the first occurrence (the entry key).
    return BIB_ENTRY_RE.sub(f"@{entry_type}{{{citekey},", bibtex, count=1)


def json_load_maybe(payload: bytes) -> dict[str, Any] | None:
    try:
        return json.loads(payload.decode("utf-8"))
    except Exception:
        return None


def first_str(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, list) and value and isinstance(value[0], str):
        return value[0]
    return ""


def main() -> int:
    parser = argparse.ArgumentParser(description="Build references.bib from a DOI list and log verification.")
    parser.add_argument("--doi-list", default="docs/references/doi_list.tsv")
    parser.add_argument("--out-bib", default="docs/references/references.bib")
    parser.add_argument("--out-verify", default="docs/CITATION_VERIFICATION.tsv")
    parser.add_argument("--cache-dir", default="docs/references/cache")
    parser.add_argument("--timeout-sec", type=float, default=20.0)
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent.parent
    doi_list_path = repo_root / args.doi_list
    out_bib_path = repo_root / args.out_bib
    out_verify_path = repo_root / args.out_verify
    cache_dir = repo_root / args.cache_dir

    rows = read_tsv(doi_list_path)
    if not rows:
        raise SystemExit(f"Empty DOI list: {doi_list_path}")

    seen: set[str] = set()
    cleaned: list[dict[str, str]] = []
    for r in rows:
        citekey = (r.get("citekey") or "").strip()
        doi = normalize_doi(r.get("doi") or "")
        notes = (r.get("notes") or "").strip()
        if not citekey or not doi:
            continue
        if citekey in seen:
            raise SystemExit(f"Duplicate citekey in DOI list: {citekey}")
        seen.add(citekey)
        if not DOI_RE.match(doi):
            raise SystemExit(f"DOI failed basic validation for {citekey}: {doi}")
        cleaned.append({"citekey": citekey, "doi": doi, "notes": notes})

    cache_dir.mkdir(parents=True, exist_ok=True)

    bib_entries: list[str] = []
    verify_rows: list[dict[str, Any]] = []
    now = utc_now_iso()

    ua = "crc-spatial-benchmark/1.0 (reference-builder)"
    def process_one(r: dict[str, str]) -> tuple[str, dict[str, Any]]:
        citekey = r["citekey"]
        doi = r["doi"]
        doi_url = f"https://doi.org/{urllib.parse.quote(doi, safe=':/')}"

        # --- BibTeX from doi.org
        bib_cache = cache_dir / f"{citekey}.bib"
        if bib_cache.exists():
            bibtex_payload = bib_cache.read_bytes()
            bib_fetch = FetchResult(
                ok=True,
                status=200,
                url=doi_url,
                final_url=doi_url,
                content_type="application/x-bibtex (cache)",
                body=bibtex_payload,
                error="",
            )
        else:
            bib_fetch = fetch_url(
                doi_url,
                headers={"Accept": "application/x-bibtex", "User-Agent": ua},
                timeout_sec=args.timeout_sec,
            )
            if bib_fetch.ok and bib_fetch.body:
                bib_cache.write_bytes(bib_fetch.body)

        bibtex_text = ""
        bibtex_ok = False
        if bib_fetch.ok and bib_fetch.body:
            try:
                bibtex_text = bib_fetch.body.decode("utf-8", errors="replace")
                bibtex_text = rewrite_bibtex_key(bibtex_text, citekey)
                bibtex_ok = bibtex_text.lstrip().startswith("@")
            except Exception:
                bibtex_ok = False

        crossref_url = f"https://api.crossref.org/works/{urllib.parse.quote(doi)}"
        openalex_url = f"https://api.openalex.org/works/https://doi.org/{urllib.parse.quote(doi)}"
        crossref_fetch = fetch_url(
            crossref_url,
            headers={"Accept": "application/json", "User-Agent": ua},
            timeout_sec=args.timeout_sec,
        )
        openalex_fetch = fetch_url(
            openalex_url,
            headers={"Accept": "application/json", "User-Agent": ua},
            timeout_sec=args.timeout_sec,
        )

        crossref_ok = False
        crossref_title = ""
        crossref_year = ""
        crossref_container = ""
        if crossref_fetch.ok and crossref_fetch.body:
            payload = json_load_maybe(crossref_fetch.body)
            msg = (payload or {}).get("message") if isinstance(payload, dict) else None
            if isinstance(msg, dict):
                crossref_ok = True
                crossref_title = first_str(msg.get("title"))
                crossref_container = first_str(msg.get("container-title"))
                issued = msg.get("issued")
                if isinstance(issued, dict):
                    parts = issued.get("date-parts")
                    if isinstance(parts, list) and parts and isinstance(parts[0], list) and parts[0]:
                        crossref_year = str(parts[0][0])

        openalex_ok = False
        openalex_title = ""
        openalex_year = ""
        openalex_venue = ""
        oa_url = ""
        fulltext_url = ""
        if openalex_fetch.ok and openalex_fetch.body:
            payload = json_load_maybe(openalex_fetch.body)
            if isinstance(payload, dict):
                openalex_ok = True
                openalex_title = (payload.get("title") or "").strip()
                openalex_year = str(payload.get("publication_year") or "").strip()
                host_venue = payload.get("host_venue")
                if isinstance(host_venue, dict):
                    openalex_venue = (host_venue.get("display_name") or "").strip()
                    fulltext_url = (host_venue.get("url") or "").strip()
                open_access = payload.get("open_access")
                if isinstance(open_access, dict):
                    oa_url = (open_access.get("oa_url") or "").strip()
                primary_loc = payload.get("primary_location")
                if isinstance(primary_loc, dict):
                    landing = (primary_loc.get("landing_page_url") or "").strip()
                    pdf = (primary_loc.get("pdf_url") or "").strip()
                    fulltext_url = oa_url or pdf or landing or fulltext_url

        fulltext_checked_url = oa_url or fulltext_url
        fulltext_accessible = ""
        if fulltext_checked_url:
            probe = fetch_url(
                fulltext_checked_url,
                headers={"User-Agent": ua},
                timeout_sec=min(args.timeout_sec, 15.0),
            )
            fulltext_accessible = "true" if probe.ok and (probe.status is None or probe.status < 400) else "false"

        title = crossref_title or openalex_title
        year = crossref_year or openalex_year
        venue = crossref_container or openalex_venue

        verify_row = {
            "citekey": citekey,
            "doi": doi,
            "doi_url": doi_url,
            "doi_resolved_url": bib_fetch.final_url or doi_url,
            "bibtex_ok": "true" if bibtex_ok else "false",
            "bibtex_status": bib_fetch.status if bib_fetch.status is not None else "",
            "crossref_ok": "true" if crossref_ok else "false",
            "crossref_status": crossref_fetch.status if crossref_fetch.status is not None else "",
            "openalex_ok": "true" if openalex_ok else "false",
            "openalex_status": openalex_fetch.status if openalex_fetch.status is not None else "",
            "two_source_verified": "true" if (bibtex_ok and (crossref_ok or openalex_ok)) else "false",
            "title": title,
            "year": year,
            "venue": venue,
            "oa_url": oa_url,
            "fulltext_url_checked": fulltext_checked_url,
            "fulltext_accessible": fulltext_accessible,
            "full_text_read": "pending",
            "full_text_read_notes": "",
            "verified_at_utc": now,
            "notes": r.get("notes", ""),
            "errors": "; ".join(x for x in [bib_fetch.error, crossref_fetch.error, openalex_fetch.error] if x),
        }

        bib_entry = (bibtex_text.strip() + "\n") if bibtex_ok else ""
        return bib_entry, verify_row

    bib_by_key: dict[str, str] = {}
    verify_by_key: dict[str, dict[str, Any]] = {}

    max_workers = min(12, max(4, (os.cpu_count() or 8)))
    with cf.ThreadPoolExecutor(max_workers=max_workers) as exec_all:
        futures = [exec_all.submit(process_one, r) for r in cleaned]
        for fut in cf.as_completed(futures):
            bib_entry, verify_row = fut.result()
            citekey = str(verify_row.get("citekey") or "")
            if citekey:
                verify_by_key[citekey] = verify_row
                bib_by_key[citekey] = bib_entry

    verify_rows = [verify_by_key[k] for k in sorted(verify_by_key)]
    bib_entries = [bib_by_key[k] for k in sorted(bib_by_key) if bib_by_key[k]]

    out_bib_path.parent.mkdir(parents=True, exist_ok=True)
    out_bib_path.write_text("\n".join(bib_entries).strip() + "\n", encoding="utf-8")

    verify_fields = [
        "citekey",
        "doi",
        "doi_url",
        "doi_resolved_url",
        "bibtex_ok",
        "bibtex_status",
        "crossref_ok",
        "crossref_status",
        "openalex_ok",
        "openalex_status",
        "two_source_verified",
        "title",
        "year",
        "venue",
        "oa_url",
        "fulltext_url_checked",
        "fulltext_accessible",
        "full_text_read",
        "full_text_read_notes",
        "verified_at_utc",
        "notes",
        "errors",
    ]
    write_tsv(out_verify_path, verify_rows, fieldnames=verify_fields)

    n_ok = sum(1 for r in verify_rows if r.get("bibtex_ok") == "true")
    n_total = len(verify_rows)
    if n_ok != n_total:
        print(f"WARNING: BibTeX fetch succeeded for {n_ok}/{n_total} references", file=sys.stderr)
        bad = [r["citekey"] for r in verify_rows if r.get("bibtex_ok") != "true"]
        print("BibTeX failures: " + ", ".join(bad[:20]) + (" ..." if len(bad) > 20 else ""), file=sys.stderr)
        return 2

    print(f"Wrote {n_total} references to {out_bib_path}")
    print(f"Wrote verification log to {out_verify_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

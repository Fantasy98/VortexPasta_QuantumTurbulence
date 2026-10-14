"""Check report links and freeze artifact provenance without running simulations."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit, unquote
from datetime import datetime, timezone
import hashlib
import json
import tomllib

ROOT = Path(__file__).resolve().parent

class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.paths = set()

    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if key in ('href', 'src') and value:
                url = urlsplit(value)
                if not url.scheme and url.path:
                    self.paths.add(unquote(url.path))

def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def main():
    final_path = ROOT / 'analysis/final_provenance.json'
    manifest = ROOT / 'SHA256SUMS'
    assert not final_path.exists() and not manifest.exists(), 'Preserve prior freeze records'
    checks = json.loads((ROOT / 'analysis/validation_checks.json').read_text())
    runs = sorted((ROOT / 'runs').glob('*/complete.toml'))
    assert len(runs) == len(checks['runs']) == 34
    for path in runs:
        done = tomllib.loads(path.read_text())
        assert done['status'] == 'complete'
        record = checks['runs'][path.parent.name]
        assert all(record[k] for k in ('finite', 'event_count_consistent',
                                      'checkpoint_retained', 'target_reached'))
    parser = Links()
    parser.feed((ROOT / 'REPORT.html').read_text())
    generated = {'analysis/final_provenance.json', 'SHA256SUMS'}
    missing = sorted(p for p in parser.paths if p not in generated and not (ROOT / p).is_file())
    assert not missing, missing
    provenance = json.loads((ROOT / 'analysis/provenance.json').read_text())
    source_paths = sorted(set(ROOT.glob('*.py')) | set(ROOT.glob('*.jl')) |
                          set(ROOT.glob('*.sh')) | set((ROOT / 'sources').glob('*.jl')) |
                          {ROOT / 'Project.toml', ROOT / 'Manifest.toml',
                           ROOT / 'requirements-analysis.txt', ROOT / 'REPRODUCE.txt'})
    metadata = {
        'frozen_utc': datetime.now(timezone.utc).isoformat(),
        'completed_runs': len(runs),
        'report_local_links_checked': len(parser.paths),
        'missing_report_local_links': missing,
        'julia': provenance['julia'],
        'VortexPasta': provenance['VortexPasta'],
        'source_groups': provenance['source_groups'],
        'source_copy_note': provenance['source_copy_note'],
        'source_sha256': {str(p.relative_to(ROOT)): digest(p) for p in source_paths},
        'manifest': 'SHA256SUMS',
        'manifest_scope': 'All regular files except SHA256SUMS and Python bytecode caches at freeze time. '
                          'The earlier analysis/provenance.json is retained unchanged; this record hashes final sources.'
    }
    final_path.write_text(json.dumps(metadata, indent=2) + '\n')
    paths = sorted(p for p in ROOT.rglob('*') if p.is_file() and p != manifest
                   and '__pycache__' not in p.parts)
    manifest.write_text(''.join(f'{digest(p)}  {p.relative_to(ROOT)}\n' for p in paths))
    print(json.dumps({'completed_runs': len(runs), 'local_links_checked': len(parser.paths),
                      'checksummed_files': len(paths),
                      'storage_MiB': sum(p.stat().st_size for p in ROOT.rglob('*') if p.is_file()) / 2**20}, indent=2))

if __name__ == '__main__':
    main()

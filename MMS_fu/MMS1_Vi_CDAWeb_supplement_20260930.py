"""取回本次核查发现的CDAWeb缺失CDF，保留原文件名并归档到Z盘；不重画图。"""
import datetime as dt
import hashlib
import importlib.util
import json
import re
import time
import urllib.parse
from pathlib import Path

code = Path(__file__).parent
spec = importlib.util.spec_from_file_location('existing_mms_download', code / 'MMS_event_overview_20260923_download.py')
d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d)
audit = d.ROOT / 'derived' / 'MMS1_tail_20260720_0815' / 'Vi_completeness_20260930'
source = json.loads((audit / 'cdaweb' / 'summary.json').read_text('utf-8'))
files = source['missing_requested_files']
assert len(files) == 49, '应仅处理本次核查发现的49个文件'
records = []
for index, item in enumerate(files, 1):
    name = item['filename']
    assert re.fullmatch(r'mms1_fpi_fast_l2_dis-moms_202608\d{8}_v3\.4\.0\.cdf', name), name
    url = item['source_url']
    parsed = urllib.parse.urlparse(url)
    assert parsed.scheme == 'https' and parsed.netloc == 'cdaweb.gsfc.nasa.gov'
    assert parsed.path.rsplit('/', 1)[-1] == name and not parsed.query
    target = d.local_path(name)
    assert target == Path(item['local_path'])
    expected = int(item['website_length'])
    if target.exists():
        assert target.stat().st_size == expected, f'已有文件大小不同，保留原文件：{target}'
        status = 'already_present'
    else:
        assert target.parent.is_dir(), target.parent
        temporary = target.with_suffix('.cdf.part')
        for attempt in range(4):
            try:
                with d.open_url(url) as response, temporary.open('wb') as output:
                    headers = dict(response.headers)
                    if response.headers.get('Content-Length'):
                        assert int(response.headers['Content-Length']) == expected, (name, 'HTTP length differs from catalogue')
                    while chunk := response.read(1024 * 1024):
                        output.write(chunk)
                assert temporary.stat().st_size == expected, (name, temporary.stat().st_size, expected)
                with temporary.open('rb') as stream:
                    assert stream.read(4) in (bytes.fromhex('cdf30001'), bytes.fromhex('cdf26002'), bytes.fromhex('0000ffff')), name
                break
            except Exception as exc:
                if attempt == 3:
                    raise
                print('RETRY', name, repr(exc), flush=True)
                time.sleep(2 * (attempt + 1))
        temporary.rename(target)
        status = 'downloaded'
    record = {**item, 'path': str(target), 'download_status': status,
              'downloaded_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
              'bytes_verified': expected, 'sha256': hashlib.sha256(target.read_bytes()).hexdigest()}
    records.append(record)
    d.save(audit / 'cdaweb_download_verification.json', {'complete': len(records) == len(files), 'files': records})
    print(f'CDF {index}/{len(files)} {status} {name}', flush=True)
print('CDAWEB_FILES_VERIFIED', len(records), sum(r['bytes_verified'] for r in records), flush=True)

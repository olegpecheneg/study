from pathlib import Path
from typing import Dict


def read_pipeline_config() -> Dict[str, str]:
    """Read simple key=value pipeline_config.conf placed in project `config/`.

    Returns dict of keys to values (quotes stripped). Non-fatal on errors.
    """
    cfg = {}
    try:
        config_path = Path(__file__).resolve().parents[2] / 'config' / 'pipeline_config.conf'
        if config_path.exists():
            text = config_path.read_text(encoding='utf-8')
            for line in text.splitlines():
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if '=' in line:
                    k, v = line.split('=', 1)
                    k = k.strip()
                    v = v.strip().strip('"').strip("'")
                    cfg[k] = v
    except Exception:
        # non-fatal: return empty config
        pass
    return cfg


def resolve_project_paths() -> Dict[str, str]:
    """Resolve common project paths with precedence: config.WORK_DIR -> project-relative defaults."""
    cfg = read_pipeline_config()
    project_root = None
    if 'WORK_DIR' in cfg and cfg['WORK_DIR']:
        try:
            project_root = Path(cfg['WORK_DIR']).expanduser().resolve()
        except Exception:
            project_root = None
    if project_root is None:
        project_root = Path(__file__).resolve().parents[2]

    paths: Dict[str, str] = {}
    paths['PROJECT_ROOT'] = str(project_root)
    paths['CONSTRUCTS_DIR'] = cfg.get('CONSTRUCTS_DIR', str(project_root / 'MitoFragilityScore' / 'Constructs'))
    paths['ENERGIES_DIR'] = cfg.get('ENERGIES_DIR', str(project_root / 'MitoFragilityScore' / 'Energies'))
    paths['OUTPUT_DIR'] = cfg.get('OUTPUT_DIR', str(project_root / 'DataPreparing' / 'output'))
    paths['SEQUENCES_RELATIVE_DIR'] = cfg.get('SEQUENCES_RELATIVE_DIR', str(project_root / 'MitoFragilityScore' / 'Sequences' / 'Relative'))
    if 'SNV_CSV' in cfg:
        p = Path(cfg['SNV_CSV'])
        if not p.is_absolute():
            p = project_root / p
        paths['SNV_CSV'] = str(p)
    # numeric defaults
    try:
        paths['DEFAULT_FDR'] = float(cfg.get('DEFAULT_FDR', 0.056))
    except Exception:
        paths['DEFAULT_FDR'] = 0.056
    try:
        paths['OUTLIER_STD_MULTIPLIER'] = float(cfg.get('OUTLIER_STD_MULTIPLIER', 2.0))
    except Exception:
        paths['OUTLIER_STD_MULTIPLIER'] = 2.0
    # JOB_ID sanitize
    raw_job = cfg.get('JOB_ID', '')
    if raw_job:
        import re
        m = re.match(r'^([A-Za-z0-9_.-]+)$', raw_job)
        if m:
            paths['JOB_ID'] = m.group(1)
        else:
            paths['JOB_ID'] = ''
    else:
        paths['JOB_ID'] = ''

    return paths

"""Locate retained production data and compacted validation records."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / 'analysis' / 'validation_runs'
ARCHIVED_STAGES = ('agrl_probe', 'agrl', 'gpe_trial', 'gpe_dt_half',
                   'reconnection', 'reconnection_dt_half')


def resolve_case_path(path):
    path = Path(path)
    if not path.is_absolute():
        path = HERE / path
    try:
        relative = path.relative_to(HERE)
    except ValueError:
        return path
    if relative.parts and relative.parts[0] in ARCHIVED_STAGES and not path.exists():
        return ARCHIVE / relative
    return path


def stage_path(name):
    return resolve_case_path(HERE / name)


def is_compacted(stage):
    return (resolve_case_path(stage) / 'cleanup_state.json').exists()

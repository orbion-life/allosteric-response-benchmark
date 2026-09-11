"""Optional, correctly scaled process metadata; never used in scientific results."""
import sys

def usage():
    try:
        import resource
        r = resource.getrusage(resource.RUSAGE_SELF)
        return {"max_rss_bytes": int(r.ru_maxrss * (1 if sys.platform == "darwin" else 1024)),
                "user_cpu_seconds": r.ru_utime, "system_cpu_seconds": r.ru_stime}
    except (ImportError, AttributeError):
        return {"max_rss_bytes": None, "user_cpu_seconds": None, "system_cpu_seconds": None}

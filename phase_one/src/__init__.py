"""
src/__init__.py — StepWise-Q
Qiskit Fall Fest 2026 | SRM University-AP | PS-S3
"""

__version__ = "0.1.0"
__all__ = ["problem", "processors", "routing", "metrics", "visualization", "recommend_steps"]


def __getattr__(name: str):
    if name == "recommend_steps":
        from src.metrics import recommend_steps
        return recommend_steps
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

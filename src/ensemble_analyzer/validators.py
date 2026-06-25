import re


def validate_line(text: str, pattern: str, expected: float, threshold: float) -> bool:
    m = re.search(pattern, text)
    if not m:
        return False
    try:
        value = float(m.group(1))
    except (ValueError, IndexError):
        return False
    return abs(value - expected) <= threshold

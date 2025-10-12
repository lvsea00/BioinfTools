def check_gc_count(seq: str, gc_lower: float, gc_upper: float) -> bool:
    """
    Checks if GC percentage of the read matches given conditions
    """
    gc_content = (seq.count('G') + seq.count('C')) * 100 / len(seq)
    return gc_lower <= gc_content <= gc_upper


def mean_quality(seq: str) -> float:
    """Counts mean quality of the read"""
    quality_sum: int = 0
    for char in seq:
        quality_sum += ord(char) - 33
    mean_quality = quality_sum/len(seq)
    return mean_quality

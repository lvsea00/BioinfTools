def check_gc_count(seq: str, gc_lower: float, gc_upper: float) -> bool:
    """
    Checks if GC percentage of the read matches given conditions
    """
    gc_content = (seq.count('G') + seq.count('C')) * 100 / len(seq)
    return gc_lower <= gc_content <= gc_upper


def check_quality(seq: str, quality_threshold: int) -> bool:
    """
    Checks if quality of the read matches given conditions
    """
    quality_sum: int = 0
    for char in seq:
        quality_sum += ord(char) - 33
    mean_quality = quality_sum/len(seq)
    return mean_quality > quality_threshold

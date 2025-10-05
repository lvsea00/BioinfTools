def count_gc(seq: str) -> float:
    """Counts GC percentage of the read"""
    gc_content = (seq.count('G') + seq.count('C')) * 100 / len(seq)
    return gc_content


def mean_quality(seq: str) -> float:
    """Counts mean quality of the read"""
    quality_sum: int = 0
    for char in seq:
        quality_sum += ord(char) - 33
    mean_quality = quality_sum/len(seq)
    return mean_quality
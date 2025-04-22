bedtools intersect -a both_with_calls_regions.bed -b GRCh38_alllowmapandsegdupregions.bed -u > intersect_hard_to_map.bed


bedtools intersect -a ont_agree_regions.bed -b GRCh38_alllowmapandsegdupregions.bed -u > intersect_ont_agree.bed

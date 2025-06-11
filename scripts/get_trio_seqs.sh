# Trio 1

#./scripts/extract_region.sh "h2tg000003l:7723021-7733021_+" data/trio1/child/chr15-20000000-50000000_child_h2.fa > breakpoint_lab/trio1/child_bp.fa
#./scripts/extract_region.sh "utg000179l:1121501-1131500_-" data/trio1/father/chr15-20000000-50000000_father_putigs.fa > breakpoint_lab/trio1/father_bp1.fa
#./scripts/extract_region.sh "utg000181l:861501-871500_-" data/trio1/father/chr15-20000000-50000000_father_putigs.fa > breakpoint_lab/trio1/father_bp2.fa

#cat breakpoint_lab/trio1/child_bp.fa breakpoint_lab/trio1/father_bp1.fa breakpoint_lab/trio1/father_bp2.fa > breakpoint_lab/trio1/trio1_all.fa
#sed -i '' '1s/.*/>child/' breakpoint_lab/trio1/trio1_all.fa # apparently works on mac only

# Trio 2 (experiment...)
# ./scripts/extract_region.sh "seq1:6375-107375_+" /Users/whops/projects/mini-programs/contig-merger/test.fa > breakpoint_lab/trio2/child_bp.fa
# ./scripts/extract_region.sh "utg000215l:653001-754000_-" data/trio2/mother/chr15-*.fa > breakpoint_lab/trio2/mother_bp2.fa
# ./scripts/extract_region.sh "utg000217l:361001-461256_-" data/trio2/mother/chr15-*.fa > breakpoint_lab/trio2/mother_bp1.fa
# cat breakpoint_lab/trio2/child_bp.fa breakpoint_lab/trio2/*_bp1.fa breakpoint_lab/trio2/*_bp2.fa > breakpoint_lab/trio2/trio2_all.fa
# sed -i '' '1s/.*/>child/' breakpoint_lab/trio2/trio2_all.fa 

# Trio 3 (20kbp)
# ./scripts/extract_region.sh "h2tg000003l:8136459-8156459_+" data/trio3/child/chr15-*.fa > breakpoint_lab/trio3/child_bp.fa
# ./scripts/extract_region.sh "utg000271l:46501-66500_+" data/trio3/father/chr15-*.fa > breakpoint_lab/trio3/father_bp1.fa
# ./scripts/extract_region.sh "utg000272l:260501-280500_-" data/trio3/father/chr15-*.fa > breakpoint_lab/trio3/father_bp2.fa
# cat breakpoint_lab/trio3/child_bp.fa breakpoint_lab/trio3/*_bp1.fa breakpoint_lab/trio3/*_bp2.fa > breakpoint_lab/trio3/trio3_all.fa
# sed -i '' '1s/.*/>child/' breakpoint_lab/trio3/trio3_all.fa 

# Trio 3 (50kbp)
# ./scripts/extract_region.sh "h2tg000003l:8121459-8171459_+" data/trio3/child/chr15-*.fa > breakpoint_lab/trio3_50/child_bp.fa
# ./scripts/extract_region.sh "utg000271l:31501-81500_+" data/trio3/father/chr15-*.fa > breakpoint_lab/trio3_50/father_bp1.fa
# ./scripts/extract_region.sh "utg000272l:245501-295500_-" data/trio3/father/chr15-*.fa > breakpoint_lab/trio3_50/father_bp2.fa
# cat breakpoint_lab/trio3_50/child_bp.fa breakpoint_lab/trio3_50/*_bp1.fa breakpoint_lab/trio3_50/*_bp2.fa > breakpoint_lab/trio3_50/trio3_all.fa
# sed -i '' '1s/.*/>child/' breakpoint_lab/trio3_50/trio3_all.fa 

# Trio 4 (3kbp)
# ./scripts/extract_region.sh "h1tg000014l:105384-108384_+" data/trio4/child/chr15-*.fa > breakpoint_lab/trio4/child_bp.fa
# ./scripts/extract_region.sh "utg000375l:24001-27000_-" data/trio4/mother/chr15-*.fa > breakpoint_lab/trio4/father_bp1.fa
# ./scripts/extract_region.sh "utg000355l:1-3000_-" data/trio4/mother/chr15-*.fa > breakpoint_lab/trio4/father_bp2.fa
# cat breakpoint_lab/trio4/child_bp.fa breakpoint_lab/trio4/*_bp1.fa breakpoint_lab/trio4/*_bp2.fa > breakpoint_lab/trio4/trio4_all.fa
# sed -i '' '1s/.*/>child/' breakpoint_lab/trio4/trio4_all.fa 

# Trio 6 (50kbp)
./scripts/extract_region.sh "h1tg000004l:7998317-8048317_+" data/trio6/child/chr15-*.fa > breakpoint_lab/trio6/child_bp.fa
./scripts/extract_region.sh "utg000348l:501501-551500_+" data/trio6/mother/chr15-*.fa > breakpoint_lab/trio6/father_bp1.fa
./scripts/extract_region.sh "utg000354l:70501-120500_+" data/trio6/mother/chr15-*.fa > breakpoint_lab/trio6/father_bp2.fa
cat breakpoint_lab/trio6/child_bp.fa breakpoint_lab/trio6/*_bp1.fa breakpoint_lab/trio6/*_bp2.fa > breakpoint_lab/trio6/trio6_all.fa
sed -i '' '1s/.*/>child/' breakpoint_lab/trio6/trio6_all.fa 
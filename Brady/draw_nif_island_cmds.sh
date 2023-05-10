
# setting
#mkdir /home-user/wzsong/Scripts/searchInBatch
#cd /home-user/wzsong/Scripts/searchInBatch
#cp -r /home-user/sswang/project/Rhizobiales/scripts/searchBlast/lib/ ./
#cp -r /home-user/sswang/project/Rhizobiales/scripts/searchBlast/lib/searchInBatch.rb ./


# blast against KEGG database using DIAMOND 
#conda activate ruby
#export RUBYLIB=/home-user/sswang/program/ruby/lib/ruby/gems/3.0.0/gems/triez-1.0.6/lib/:$RUBYLIB
cd /home-user/wzsong/draw_nif_island
ruby /home-user/wzsong/Scripts/searchInBatch/searchInBatch.rb --indir protein/ --outdir kegg_diamond/ --cpu 4 --kegg --diamond --force --thread 4


# get seq_info.tbl
cd /home-user/wzsong/draw_nif_island
ruby /home-user/sswang/project/Rhizobiales/scripts/correlTraitGene/parseSearchRes.rb -e 1e-10 --indir kegg_diamond/ --outdir kegg_res/ --force --db_type kegg --cpu 4 --fam_list /home-user/sswang/resource/db/kegg/parse/koid_all.2017.tbl --seq_indir protein/ --output_seq


# add header
cd /home-user/wzsong/draw_nif_island
ruby /home-user/sswang/project/Brady/scripts/geneArrangement_Jinjin/get_header.rb -i seq_info.tbl > seq_info_header.tbl


# add the header line in seq_info_header.txt to seq_info.tbl
cat seq_info_header.tbl seq_info.tbl > seq_info_with_header.tbl


# Tao's script:
cd /home-user/wzsong/draw_nif_island
python3 /home-user/sswang/project/Brady/scripts/geneArrangement_Jinjin/GeneArrangementFigure.py --gbk_indir gbk/ --annotation seq_info_with_header.tbl --species_list species_list.txt --length 70000
python3 /home-user/sswang/project/Brady/scripts/geneArrangement_Jinjin/GeneArrangementFigure.py --gbk_indir gbk/ --annotation seq_info_with_header.tbl --species_list species_list.txt --length 70000 --one





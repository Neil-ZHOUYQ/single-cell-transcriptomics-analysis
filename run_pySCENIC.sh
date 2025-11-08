pip install pyscenic==0.12.1
pip install numpy==1.23.5
pip install pandas==1.5.3
pip install numba==0.56.4
pip install dask==2024.5.0
pip install dask-expr==1.1.0
conda install -y scanpy==1.10.0


# #查看所有镜像
# conda config --show channels

# #删除全部镜像
# conda config --remove-key channels

# #添加镜像
# conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
# conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
# conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
# conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
# conda config --add channels bioconda
# conda config --add channels conda-forge
# conda config --add channels https://mirrors.ustc.edu.cn/anaconda/pkgs/main/
# conda config --add channels https://mirrors.ustc.edu.cn/anaconda/pkgs/free/
# conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/conda-forge/
# conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/msys2/
# conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/bioconda/
# conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/menpo/
# conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/
# conda config --add channels r 
# conda config --add channels conda-forge 
# conda config --add channels bioconda
# conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda/
# conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge/
# conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/free/
# conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/main/
# conda config --set show_channel_urls yes
# conda config --set always_yes True

# #创建环境
# conda create -y -n pyscenic_env python=3.10

# #激活环境
# conda activate pyscenic_env



#人
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
wget https://github.com/aertslab/pySCENIC/blob/master/resources/hs_hgnc_tfs.txt

#鼠
wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl
wget https://github.com/aertslab/pySCENIC/blob/master/resources/mm_mgi_tfs.txt





python create_loom.py








pyscenic grn \
--num_workers 40 \
--sparse \
--method grnboost2 \
--output sce.adj.csv \
sce.loom \
hs_hgnc_tfs.txt

pyscenic ctx \
--num_workers 40 \
--output sce.regulons.csv \
--expression_mtx_fname sce.loom \
--all_modules \
--mask_dropouts \
--mode "dask_multiprocessing" \
--min_genes 10 \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
sce.adj.csv \
hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather

pyscenic aucell \
--num_workers 40 \
--output sce_SCENIC.loom \
sce.loom \
sce.regulons.csv
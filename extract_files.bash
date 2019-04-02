mkdir acc blca brca cesc chol coadread dlbc esca gbm hnsc kich kirc kirp laml lgg lihc luad lusc meso ov paad pcpg prad sarc skcm stad tgct thca thym ucec ucs uvm

#ls *tar.gz|perl -pe 's/(.+)_tcga_pan.+/tar -xzf $1_tcga_pan_can_atlas_2018.tar.gz -C $1\//g'

tar -xzf ov_tcga_pan_can_atlas_2018.tar.gz -C ov/

tar -xzf acc_tcga_pan_can_atlas_2018.tar.gz -C acc/
tar -xzf blca_tcga_pan_can_atlas_2018.tar.gz -C blca/
tar -xzf brca_tcga_pan_can_atlas_2018.tar.gz -C brca/
tar -xzf cesc_tcga_pan_can_atlas_2018.tar.gz -C cesc/
tar -xzf chol_tcga_pan_can_atlas_2018.tar.gz -C chol/
tar -xzf coadread_tcga_pan_can_atlas_2018.tar.gz -C coadread/
tar -xzf dlbc_tcga_pan_can_atlas_2018.tar.gz -C dlbc/
tar -xzf esca_tcga_pan_can_atlas_2018.tar.gz -C esca/
tar -xzf gbm_tcga_pan_can_atlas_2018.tar.gz -C gbm/
tar -xzf hnsc_tcga_pan_can_atlas_2018.tar.gz -C hnsc/
tar -xzf kich_tcga_pan_can_atlas_2018.tar.gz -C kich/
tar -xzf kirc_tcga_pan_can_atlas_2018.tar.gz -C kirc/
tar -xzf kirp_tcga_pan_can_atlas_2018.tar.gz -C kirp/
tar -xzf laml_tcga_pan_can_atlas_2018.tar.gz -C laml/
tar -xzf lgg_tcga_pan_can_atlas_2018.tar.gz -C lgg/
tar -xzf lihc_tcga_pan_can_atlas_2018.tar.gz -C lihc/
tar -xzf luad_tcga_pan_can_atlas_2018.tar.gz -C luad/
tar -xzf lusc_tcga_pan_can_atlas_2018.tar.gz -C lusc/
tar -xzf meso_tcga_pan_can_atlas_2018.tar.gz -C meso/
tar -xzf paad_tcga_pan_can_atlas_2018.tar.gz -C paad/
tar -xzf pcpg_tcga_pan_can_atlas_2018.tar.gz -C pcpg/
tar -xzf prad_tcga_pan_can_atlas_2018.tar.gz -C prad/
tar -xzf sarc_tcga_pan_can_atlas_2018.tar.gz -C sarc/
tar -xzf skcm_tcga_pan_can_atlas_2018.tar.gz -C skcm/
tar -xzf stad_tcga_pan_can_atlas_2018.tar.gz -C stad/
tar -xzf tgct_tcga_pan_can_atlas_2018.tar.gz -C tgct/
tar -xzf thca_tcga_pan_can_atlas_2018.tar.gz -C thca/
tar -xzf thym_tcga_pan_can_atlas_2018.tar.gz -C thym/
tar -xzf ucec_tcga_pan_can_atlas_2018.tar.gz -C ucec/
tar -xzf ucs_tcga_pan_can_atlas_2018.tar.gz -C ucs/
tar -xzf uvm_tcga_pan_can_atlas_2018.tar.gz -C uvm/


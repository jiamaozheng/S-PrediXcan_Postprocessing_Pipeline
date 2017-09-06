#!/usr/bin/env python

from helpers import * 

__author__ = "Jiamao Zheng <jiamaoz@yahoo.com>"
__version__ = "Revision: 0.0.0.1"
__date__ = "Date: 2017-08-27"

class MetaXcanPostprocess(object):

    def __init__(self):
        # source code file path 
        self.currentPath = ''       # /src/ path 
        self.srcPath = ''

        # input, intermediate output, output file paths 
        self.output_annotation_path = '' 
        self.input_path = ''
        self.output_qqplot_path = ''
        self.manhattanplot_output_path = ''
        self.top_genes_output_path = ''
        self.top_genes_snps_output_path = ''
        self.sorted_top_genes_all_tissues_path = ''
        self.bubble_plot_output_path = ''
        self.region_plot_output_path = ''
        self.locuszoom_plot_path = ''

        # csv files 
        self.outputFileListWithoutSNPs = ''   # all metaxcan *.csv files 
        self.merged_output = ''               # merged metaxcan .csv file 

        # data 
        self.top_gene_list_dataframe = ''
        self.tissue_abbr = ''       # tissue names 
        self.gwas_lead_snp = ''     # top gene lists across tissues 
        self.data = ''              # merged metaxcan output (dataframe)
        self.total_cut_off = ''
        self.geneOfInterest = ''

        # logger 
        self.logger = logging.getLogger()

        # time 
        self.start_time = time.time() 

        # prefix for querying db 
        self.SQL_QUERY_PREFIX = "select e.genename, w.rsid from weights w join extra e on w.gene = e.gene where e.genename = '"
        self.SQL_QUERY_PREFIX_DNG = "select gene, rsid from weights where gene = '"

        # project INFO
        self.project_name = ''
        self.project_time = ''
        self.project_id = ''
        self.project = ''

        # options
        self.multiple_tissue = '' 
        self.locuszoom = ''
        self.models_folder = ''
        self.metaxcan_folder = ''
        self.tools_folder = ''


    def flushPipeline(self):
        path ='../locuszoom/locuszoom_plots/'
        os.system("rm -rf %s" % path)


    def get_args(self):
        # setup commond line arguments 
        parser = argparse.ArgumentParser()

        # required arguments
        parser.add_argument('-p', '--project_name', required=True, default=None, type=str, help='e.g breast_cancer, or multiple_tissues')

        # metaxcan outputs  
        parser.add_argument('-f', '--metaxcan_folder', required=True, default='../data/metaxcan/', type=str, help='file path to metaxcan outputs')

        # db 
        parser.add_argument('-d', '--models_folder', required=False, default='../data/models/', type=str, help='file path to prediction models')

        # plink 
        parser.add_argument('-t', '--tools_folder', required=False, default='../data/tools/', type=str, help='plink software')

        # not required, and this is optional argument. Please type 'true' if you would like to run multiple tissue pipeline 
        parser.add_argument('-m', '--multiple_tissue', required=False, default='false', type=str, help='true, if you would like to analyze outputs from multiple_tissue pipeline')

 		# not required, and this is optional argument. Please type 'false' if you don't want to run locuszoom 
        parser.add_argument('-l', '--locuszoom', required=False, default='true', type=str, help='false, if you do not want to run locuszoom')

        # parse the arguments 
        args = parser.parse_args()

        self.project_name = args.project_name 
        self.multiple_tissue = args.multiple_tissue 
        self.locuszoom = args.locuszoom
        self.models_folder = args.models_folder
        self.metaxcan_folder = args.metaxcan_folder
        self.tools_folder = args.tools_folder


    def get_parameters(self):
        self.project_time = datetime.now().strftime('%Y-%m-%d') 
        self.project_id = str(myuuid.uuid4())
        self.project = self.project_name + '_' + self.project_id + '_' + self.project_time

    # Annotations 
    def annotateMetaxcanOutput(self):
        projectName = self.project 
        # set up logger  
        log_path = '../output/' + projectName + '/logs/'
        os.makedirs(log_path)
        getLog(self.logger, log_path + projectName + '.log')

        # create directory that holds output and inputs 
        self.output_annotation_path = '../output/' + projectName  + '/annotated_output_files/' 
        os.makedirs(self.output_annotation_path)

        # self.input_path = '../input/'
        # self.input_path = '../input/' + projectName +'/'
        # os.makedirs(self.input_path)
        # os.system('mv ../input/*.* ' + self.input_path)
        # os.system('mv ../input/plink ' + self.input_path)
        self.input_path = self.metaxcan_folder


        # get input file lists (raw metaxcan output files *.csv) 
        inputFileList = glob.glob(self.input_path + '*.csv')

        # loop through file lists 
        for inputFilename in inputFileList:
            msg = "\n " + datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "ANNOTATING GENES: " + inputFilename 
            self.logger.info(msg)
            print(msg)
            # read data 
            r("data <- fread('%s')" %inputFilename)
            data = r('setDF(data)')
            data = r('na.omit(data)')
            robjects.globalenv['dataframe'] = data

            # annotation library 
            grch38 = r('grch38')
            robjects.globalenv['dataframe'] = grch38

            # annotating  
            # if inputFilename == 'xxxxxDGN-WB-unscaled.csv':
            #     annotatedData = r("inner_join(data[, c('gene', 'gene_name', 'zscore', 'effect_size', 'pvalue', 'n_snps_in_model', 'pred_perf_p', 'pred_perf_R2')], grch37[, c('symbol', 'chr', 'start', 'end')], by=c('gene'='symbol'))")
            # else:
            annotatedData = r("inner_join(data, grch38, by =c('gene'='ensgene'))")
            annotatedData.drop(['entrez'], axis=1, inplace=True, errors='ignore')
            annotatedData = annotatedData.drop_duplicates()
            annotatedData = annotatedData[(annotatedData.biotype == 'protein_coding')]
            # annotatedData.rename(columns={'symbol' : 'gene_name'}, inplace=True)
            annotatedData['gene_name'] = annotatedData['symbol']

            # ouput annotated data 
            annotatedData.to_csv(self.output_annotation_path + inputFilename.split('/')[-1][:-4] + '_annotated.csv', index=None)

        # merge all annotated metaxcan files 
        # get annotated file lists 
        self.outputFileListWithoutSNPs = glob.glob(self.output_annotation_path + "*.csv")

        dfListWithoutSNPs = []
        for outputFileName in self.outputFileListWithoutSNPs:
            msg = "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "CONCATING FILE: " + outputFileName 
            self.logger.info(msg)
            print(msg)

            r("data <- read.csv('%s')" %outputFileName)
            data = r('data <- na.omit(data)')
            robjects.globalenv['dataframe'] = data

            # outputFileName = outputFileName.split('/')[-1]

            tissues = ['multitissue', 'TW_Liver', 'TW_Brain_Cerebellar_Hemisphere', 'TW_Esophagus_Muscularis', 'TW_Skin_Not_Sun_Exposed_Suprapubic', 'TW_Brain_Caudate_basal_ganglia', 'TW_Heart_Atrial_Appendage', 'TW_Artery_Coronary', 'TW_Esophagus_Gastroesophageal_Junction', 'TW_Adipose_Subcutaneous', 'TW_Stomach', 'TW_Artery_Tibial', 'TW_Pancreas', 'TW_Prostate', 'TW_Testis', 'TW_Brain_Cerebellum', 'TW_Vagina', 'TW_Thyroid', 'TW_Colon_Sigmoid', 'TW_Cells_Transformed_fibroblasts', 'TW_Adipose_Visceral_Omentum', 'TW_Brain_Frontal_Cortex_BA9', 'TW_Spleen', 'TW_Whole_Blood', 'TW_Brain_Hippocampus', 'TW_Pituitary', 'TW_Lung', 'TW_Brain_Nucleus_accumbens_basal_ganglia', 'TW_Esophagus_Mucosa', 'TW_Nerve_Tibial', 'TW_Heart_Left_Ventricle', 'TW_Brain_Anterior_cingulate_cortex_BA24', 'TW_Ovary', 'TW_Brain_Cortex', 'TW_Adrenal_Gland', 'TW_Muscle_Skeletal', 'TW_Cells_EBV-transformed_lymphocytes', 'TW_Artery_Aorta', 'TW_Colon_Transverse', 'TW_Breast_Mammary_Tissue', 'TW_Skin_Sun_Exposed_Lower_leg', 'TW_Brain_Putamen_basal_ganglia', 'TW_Small_Intestine_Terminal_Ileum', 'TW_Uterus', 'TW_Brain_Hypothalamus']

            outputFileName = outputFileName.split('/')[-1]
            for tissue in tissues: 
                if tissue in outputFileName:
                    outputFileName = tissue 
            # if 'TW_' in outputFileName: 
            #     outputFileName = outputFileName[3:]

            # if 'CrossTissue_elasticNet' in outputFileName:
            #     outputFileName = outputFileName[:-25]
            # elif 'DGN-WB-unscaled' in outputFileName: 
            #     outputFileName = outputFileName [:-23]
            # else:
            #     if '_elasticNet' in outputFileName: 
            #         outputFileName = outputFileName[:-25]
            #     else: 
            #         outputFileName = outputFileName[:-14]

            # tissue = outputFileName.split('/')[-1][3:]
            data.insert(2, 'tissue', outputFileName)
            dfListWithoutSNPs.append(data)
            
        concatDfWithoutSNPs = pandas.concat(dfListWithoutSNPs, axis = 0)

        #output merged data 
        self.merged_output = self.output_annotation_path + 'merged_annotated.csv'
        concatDfWithoutSNPs.to_csv(self.merged_output, index = None)
        msg = "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done for Annotating metaxcan output!"
        self.logger.info(msg)
        print(msg)

    #  QQ-Plot  
    def createQQPlot(self):
        projectName = self.project 

        # create directory that holds qq plot output
        self.output_qqplot_path = '../output/' + projectName  + '/qq-plot/' 
        os.system("mkdir " + self.output_qqplot_path)

        # loop through file lists 
        for outputFileName in self.outputFileListWithoutSNPs: 
            msg = "\n " + datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "QQ-PLOT(no snps): " + outputFileName 
            self.logger.info(msg)
            print(msg)
            # get GWAS output data 
            r("data <- read.csv('%s')" %outputFileName)
            data = r('data <- na.omit(data)')
            robjects.globalenv['dataframe'] = data

            # draw qq-plot and save them to files
            tissues = ['multitissue', 'TW_Liver', 'TW_Brain_Cerebellar_Hemisphere', 'TW_Esophagus_Muscularis', 'TW_Skin_Not_Sun_Exposed_Suprapubic', 'TW_Brain_Caudate_basal_ganglia', 'TW_Heart_Atrial_Appendage', 'TW_Artery_Coronary', 'TW_Esophagus_Gastroesophageal_Junction', 'TW_Adipose_Subcutaneous', 'TW_Stomach', 'TW_Artery_Tibial', 'TW_Pancreas', 'TW_Prostate', 'TW_Testis', 'TW_Brain_Cerebellum', 'TW_Vagina', 'TW_Thyroid', 'TW_Colon_Sigmoid', 'TW_Cells_Transformed_fibroblasts', 'TW_Adipose_Visceral_Omentum', 'TW_Brain_Frontal_Cortex_BA9', 'TW_Spleen', 'TW_Whole_Blood', 'TW_Brain_Hippocampus', 'TW_Pituitary', 'TW_Lung', 'TW_Brain_Nucleus_accumbens_basal_ganglia', 'TW_Esophagus_Mucosa', 'TW_Nerve_Tibial', 'TW_Heart_Left_Ventricle', 'TW_Brain_Anterior_cingulate_cortex_BA24', 'TW_Ovary', 'TW_Brain_Cortex', 'TW_Adrenal_Gland', 'TW_Muscle_Skeletal', 'TW_Cells_EBV-transformed_lymphocytes', 'TW_Artery_Aorta', 'TW_Colon_Transverse', 'TW_Breast_Mammary_Tissue', 'TW_Skin_Sun_Exposed_Lower_leg', 'TW_Brain_Putamen_basal_ganglia', 'TW_Small_Intestine_Terminal_Ileum', 'TW_Uterus', 'TW_Brain_Hypothalamus']

            outputFileName = outputFileName.split('/')[-1]
            for tissue in tissues: 
                if tissue in outputFileName:
                    outputFileName = tissue 
            # if 'TW_' in outputFileName: 
            #     outputFileName = outputFileName[3:]

            # if 'CrossTissue_elasticNet' in outputFileName:
            #     outputFileName = outputFileName[:-25]
            # elif 'DGN-WB-unscaled' in outputFileName: 
            #     outputFileName = outputFileName [:-23]
            # else:
            #     if '_elasticNet' in outputFileName: 
            #         outputFileName = outputFileName[:-25]
            #     else: 
            #         outputFileName = outputFileName[:-14] 

            r.pdf('%s%s%s%s'%(self.output_qqplot_path, 'QQ-plot_', outputFileName,'.pdf'))
            qqman.qq(data['pvalue'],  main = outputFileName)
            r['dev.off']()
        msg = "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "QQ-PLOT(no snps): " + "all_tissues"
        self.logger.info(msg)
        print(msg)

        # set up path and read merged annotated metaxcan output files  
        r("data <- read.csv('%s')" % (self.merged_output))
        data = r('data <- na.omit(data)')
        robjects.globalenv['dataframe'] = data

        # draw qq-plot and save them to files 
        r.pdf('%s%s%s%s'%(self.output_qqplot_path,'QQ-plot_', 'all_tissues', '.pdf'))
        qqman.qq(data['pvalue'],  main = "All Tissues")
        r['dev.off']()

        msg =  "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done for QQ Plot!"
        self.logger.info(msg)
        print(msg)
    
    #  Manhattan-Plot 
    def createManhattanPlot(self):
        projectName = self.project 

        # create directory that holds manhattan plot output
        self.manhattanplot_output_path = '../output/' + projectName  + '/manhattan-plot/'
        os.system("mkdir " + self.manhattanplot_output_path)

        # loop through file lists
        for outputFileName in self.outputFileListWithoutSNPs: 
            msg = "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "MANHATTAN-PLOT(no snps): " + outputFileName
            self.logger.info(msg)
            print(msg)
            # read data 
            r("data <- read.csv('%s')" %outputFileName)
            r('data$chr <- as.numeric(as.character(data$chr))')
            data = r('data <- na.omit(data)')
            robjects.globalenv['dataframe'] = data

            outputFileName = outputFileName.split('/')[-1]

            tissues = ['multitissue', 'TW_Liver', 'TW_Brain_Cerebellar_Hemisphere', 'TW_Esophagus_Muscularis', 'TW_Skin_Not_Sun_Exposed_Suprapubic', 'TW_Brain_Caudate_basal_ganglia', 'TW_Heart_Atrial_Appendage', 'TW_Artery_Coronary', 'TW_Esophagus_Gastroesophageal_Junction', 'TW_Adipose_Subcutaneous', 'TW_Stomach', 'TW_Artery_Tibial', 'TW_Pancreas', 'TW_Prostate', 'TW_Testis', 'TW_Brain_Cerebellum', 'TW_Vagina', 'TW_Thyroid', 'TW_Colon_Sigmoid', 'TW_Cells_Transformed_fibroblasts', 'TW_Adipose_Visceral_Omentum', 'TW_Brain_Frontal_Cortex_BA9', 'TW_Spleen', 'TW_Whole_Blood', 'TW_Brain_Hippocampus', 'TW_Pituitary', 'TW_Lung', 'TW_Brain_Nucleus_accumbens_basal_ganglia', 'TW_Esophagus_Mucosa', 'TW_Nerve_Tibial', 'TW_Heart_Left_Ventricle', 'TW_Brain_Anterior_cingulate_cortex_BA24', 'TW_Ovary', 'TW_Brain_Cortex', 'TW_Adrenal_Gland', 'TW_Muscle_Skeletal', 'TW_Cells_EBV-transformed_lymphocytes', 'TW_Artery_Aorta', 'TW_Colon_Transverse', 'TW_Breast_Mammary_Tissue', 'TW_Skin_Sun_Exposed_Lower_leg', 'TW_Brain_Putamen_basal_ganglia', 'TW_Small_Intestine_Terminal_Ileum', 'TW_Uterus', 'TW_Brain_Hypothalamus']

            outputFileName = outputFileName.split('/')[-1]
            for tissue in tissues: 
                if tissue in outputFileName:
                    outputFileName = tissue 
            # if 'TW_' in outputFileName: 
            #     outputFileName = outputFileName[3:]

            # if 'CrossTissue_elasticNet' in outputFileName:
            #     outputFileName = outputFileName[:-25]
            # elif 'DGN-WB-unscaled' in outputFileName: 
            #     outputFileName = outputFileName [:-23]
            # else:
            #     if '_elasticNet' in outputFileName: 
            #         outputFileName = outputFileName[:-25]
            #     else: 
            #         outputFileName = outputFileName[:-14] 

            # draw manhattan and save them to files 
            # suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
            r.pdf('%s%s%s%s'%(self.manhattanplot_output_path, 'Manhattan-plot_', outputFileName, '.pdf'))
            qqman.manhattan(data, chr = 'chr', bp='start', p='pvalue', snp='gene_name', cex = 0.5, suggestiveline = 'FALSE', genomewideline = 'FALSE', main = outputFileName)
            r['dev.off']()

        # read data 
        r("data <- read.csv('%s')" % self.merged_output)
        r('data$chr <- as.numeric(as.character(data$chr))')
        data = r('data <- na.omit(data)')
        robjects.globalenv['dataframe'] = data

        r("snpOfInterest <- read.csv('%s')" % self.sorted_top_genes_all_tissues_path)
        self.geneOfInterest = r('snpOfInterest <- as.vector(snpOfInterest$gene_name)')
        robjects.globalenv['vector'] = self.geneOfInterest

        # print(self.geneOfInterest)

        # draw manhattan and save them to files 
        r.pdf('%s%s%s%s'%(self.manhattanplot_output_path, 'manhattan-plot_', 'all_tissues', '.pdf'))
        # snpOfInterest = self.top_gene_list_dataframe[['gene_name']]
        # print(snpOfInterest)
        # highlight = self.geneOfInterest,
        qqman.manhattan(data, chr = 'chr', bp='start', p='pvalue', snp='gene_name', cex = 0.5, suggestiveline = 'FALSE', genomewideline = 'FALSE', annotatePval = self.total_cut_off, annotateTop = 'FALSE', main = 'All Tissues')
        r['dev.off']()

        msg =  "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done for Manhattan Plot!"
        self.logger.info(msg)
        print(msg)

    # Top Gene List without SNPs
    def getTopGeneList(self):
        projectName = self.project 

        # create directory that holds top gene list (no associated snp)
        self.top_genes_output_path = '../output/' + projectName  +'/top_genes/'
        os.system("mkdir " + self.top_genes_output_path)

        # loop through files 
        for outputFileName in self.outputFileListWithoutSNPs:
            msg = "\n " + datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "TOP GENE LIST: " + outputFileName
            self.logger.info(msg)
            print(msg)

            df = pandas.read_csv(outputFileName) 

            tissues = ['multitissue', 'TW_Liver', 'TW_Brain_Cerebellar_Hemisphere', 'TW_Esophagus_Muscularis', 'TW_Skin_Not_Sun_Exposed_Suprapubic', 'TW_Brain_Caudate_basal_ganglia', 'TW_Heart_Atrial_Appendage', 'TW_Artery_Coronary', 'TW_Esophagus_Gastroesophageal_Junction', 'TW_Adipose_Subcutaneous', 'TW_Stomach', 'TW_Artery_Tibial', 'TW_Pancreas', 'TW_Prostate', 'TW_Testis', 'TW_Brain_Cerebellum', 'TW_Vagina', 'TW_Thyroid', 'TW_Colon_Sigmoid', 'TW_Cells_Transformed_fibroblasts', 'TW_Adipose_Visceral_Omentum', 'TW_Brain_Frontal_Cortex_BA9', 'TW_Spleen', 'TW_Whole_Blood', 'TW_Brain_Hippocampus', 'TW_Pituitary', 'TW_Lung', 'TW_Brain_Nucleus_accumbens_basal_ganglia', 'TW_Esophagus_Mucosa', 'TW_Nerve_Tibial', 'TW_Heart_Left_Ventricle', 'TW_Brain_Anterior_cingulate_cortex_BA24', 'TW_Ovary', 'TW_Brain_Cortex', 'TW_Adrenal_Gland', 'TW_Muscle_Skeletal', 'TW_Cells_EBV-transformed_lymphocytes', 'TW_Artery_Aorta', 'TW_Colon_Transverse', 'TW_Breast_Mammary_Tissue', 'TW_Skin_Sun_Exposed_Lower_leg', 'TW_Brain_Putamen_basal_ganglia', 'TW_Small_Intestine_Terminal_Ileum', 'TW_Uterus', 'TW_Brain_Hypothalamus']

            outputFileName = outputFileName.split('/')[-1]
            for tissue in tissues: 
                if tissue in outputFileName:
                    outputFileName = tissue 
                # if 'TW_' in outputFileName: 
                #     outputFileName = outputFileName[3:]

                # if 'CrossTissue_elasticNet' in outputFileName:
                #     outputFileName = outputFileName[:-25]
                # elif 'DGN-WB-unscaled' in outputFileName: 
                #     outputFileName = outputFileName [:-23]
                # else:
                #     if '_elasticNet' in outputFileName: 
                #         outputFileName = outputFileName[:-25]
                #     else: 
                #         outputFileName = outputFileName[:-14] 

            df.insert(2, 'tissue', outputFileName)

            # sort data by defined column 
            df.sort_values(['pvalue'], ascending=[1], inplace=True)
            total_rows = df.shape[0] # number of row count shape[1] is number of col count 

            # total_rows_in_total += total_rows 
            cut_off = 0.05/total_rows 

            # cut_off_list.append(cut_off)
            top_gene_list = df[df['pvalue'] < cut_off]
            # top_gene_list = top_gene_list[top_gene_list['pred_perf_R2'] > 0.01]
            # top_gene_list.drop('gene', axis=1, inplace=True)
            msg = "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + 'CALCULATING TISSUE-WIDE P-VALUE: ' + str(cut_off)
            self.logger.info(msg)
            print(msg)

            #output data 
            top_gene_list.to_csv(self.top_genes_output_path + "sorted_top_genes_%s"%outputFileName + ".csv", index = None)

            # dfList.append(top_gene_list)

        df = pandas.read_csv(self.merged_output) 
        total_rows_in_total = df.shape[0] 
        msg = "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + 'total_row_in_total: ' + str(total_rows_in_total)
        self.logger.info(msg)
        print(msg)

        self.total_cut_off = 0.05/total_rows_in_total 
        msg =  "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + 'CALCULATING GENOME-WIDE P-VALUE: ' + str(self.total_cut_off)
        self.logger.info(msg)
        print(msg)

        df = df[df['pvalue'] < self.total_cut_off]

        # sort data by defined column 
        df.sort_values(['pvalue'], ascending=[1], inplace=True)

        #output data 
        self.sorted_top_genes_all_tissues_path = self.top_genes_output_path+ 'sorted_top_genes_all_tissues.csv' 
        df.to_csv(self.sorted_top_genes_all_tissues_path, index = None)
        self.top_gene_list_dataframe = df

        msg = "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done for Top Gene List!"
        self.logger.info(msg)
        print(msg)

    #  Top Gene List with SNPs 
    def getTopGeneListWithSNPs(self):
        projectName = self.project 

        # create directory that holds top gene list (no associated snp)
        self.top_genes_snps_output_path = '../output/' + projectName  + '/top_genes_snps/'
        os.system("mkdir " + self.top_genes_snps_output_path)

        # read data 
        top_genes = pandas.read_csv(self.sorted_top_genes_all_tissues_path)

        # add into files 
        gene_lists = top_genes['gene_name']
        tissue_lists = top_genes['tissue']
        pvalue_lists = top_genes['pvalue']
        zscore_lists = top_genes['zscore']
        model_n_lists = top_genes['n_snps_in_model']
        pred_perf_R2_lists = top_genes['pred_perf_r2']
        chr_lists = top_genes['chr']
        start_lists = top_genes['start']
        end_lists = top_genes['end']

        # get db lists 
        self.input_path = self.models_folder
        dbFileList = glob.glob(self.input_path + "*.db")

        database_names = []
        for dbFilename in dbFileList:
           database_names.append(dbFilename)

        # Loop through databases 
        query_output_list = []
        for i in range(len(database_names)):
            for k in range(len(tissue_lists)): 
                if tissue_lists[k] in database_names[i]: 
                    # print(tissue_lists[k])
                    # Connect databases 
                    conn = sqlite3.connect(database_names[i]) 

                    full_query_name = None 
                    # print(gene_lists[k])
                    if database_names[i].split('/')[-1] == 'DGN-WB-unscaled_0.5.db':
                        full_query_name = self.SQL_QUERY_PREFIX_DNG + gene_lists[k] + "'"
                    else: 
                        full_query_name = self.SQL_QUERY_PREFIX + gene_lists[k] + "'"

                    query_output = pandas.read_sql(full_query_name, conn, index_col=None)

                    if database_names[i].split('/')[-1] == 'DGN-WB-unscaled_0.5.db':
                        query_output.rename(columns={'gene':'genename'}, inplace=True) 

                    # # Add correspinding parameters to the new output file 
                    query_output['tissue'] = tissue_lists[k]
                    query_output['pvalue'] =  pvalue_lists[k]
                    query_output['zscore'] =  zscore_lists[k]
                    query_output['model_n'] = model_n_lists[k]
                    query_output['pred_perf_r2'] = pred_perf_R2_lists[k]
                    query_output['chr'] =  chr_lists[k]
                    query_output['start'] = start_lists[k]
                    query_output['end'] = end_lists[k]
                    query_output_list.append(query_output) 

                    msg = "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + 'FETCH SNPs for GENE %s ' % gene_lists[k] +  'FROM DATABASE: %s' % database_names[i]
                    self.logger.info(msg)
                    print(msg)

                    # Close database
                    conn.close() 

        # Merge output data 
        query_output_of = pandas.concat(query_output_list, axis = 0)   
        # gwas_snp = pandas.read_csv('gwas_snp.txt', sep="\s+")
        # gwas_snp.rename(columns={'MarkerName':'rsid', 'P-value':'gwas_pvalue'}, inplace=True)

        # query_output_of.merge(gwas_snp, on='rsid', how='inner')

        # Output merged data 
        query_output_of.to_csv(self.top_genes_snps_output_path + "top_genes_snps.csv", index=None)

        msg = "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done for Top Gene List with SNPs!"
        self.logger.info(msg)
        print(msg)

    #  Bubble Plot 
    def createBubblePlot(self):
        projectName = self.project 

        # create directory that holds top gene list (no associated snp)
        self.bubble_plot_output_path = '../output/' + projectName  + '/bubble_plot/'
        os.system("mkdir " + self.bubble_plot_output_path)

        # current path 
        self.srcPath = os.getcwd()
        self.currentPath = os.getcwd()[:-((len(os.getcwd().split('/')[-1])) + 1)]
        # print(os.getcwd().split('/')[-1])
        # print(self.currentPath)
        # print(self.currentPath + '/input/' + projectName +'/')
        # os.chdir(self.currentPath + '/input/' + projectName +'/')

        # get tissue abbr name 
        # self.tissue_abbr = r("""

        # tissue_abbr <- read.delim('gtex_tissue_abbr.txt', sep='\t') %>% 
        # select(tissue_site_detail_abbr,tissue_site_detail_id)
        # colnames(tissue_abbr)[colnames(tissue_abbr) == 'tissue_site_detail_id'] = "tissue"

        # """) 
        # robjects.globalenv['dataframe'] = self.tissue_abbr

        # get gene lists 
        os.chdir(self.currentPath + '/output/' + projectName + '/top_genes/')
        self.gwas_lead_snp = r("""
              gwas_lead_snp <- read.csv('sorted_top_genes_all_tissues.csv') %>% 
              select(gene_name, chr, start) %>%
              distinct()
            colnames(gwas_lead_snp) = c('snpsNames', 'chrosome', 'startSites')

            """)
        robjects.globalenv['dataframe'] = self.gwas_lead_snp

        # get merged data 
        os.chdir(self.currentPath + '/output/' + projectName + '/annotated_output_files/')
        self.data = r("data <- read.csv('merged_annotated.csv')")
        robjects.globalenv['dataframe'] = self.data 

        os.chdir(self.currentPath + '/output/' + projectName + '/bubble_plot/')

        # for loop through each snps site +/- 1000000 bp  
        r("""
            # read dataframe 
            startSites = as.numeric(as.character(gwas_lead_snp$startSites))
            snpsNames = as.character(gwas_lead_snp$snpsNames)
            chrosome = as.numeric(as.character(gwas_lead_snp$chrosome))

            # for loop through each snps site +/- 1000000 bp  
            for (i in 1:length(startSites)) 
            {
                print(paste("BUBBLE PLOTING", ": ", snpsNames[i]))

                # subset data for each snps 
                subData <- subset(data, data$start > startSites[i] - 1000000  & 
                    data$start < startSites[i] + 1000000 & data$chr==chrosome[i])
                subData <- subData[order(subData$start),]
                subData <- mutate(subData, z_score=ifelse(subData$zscore > 0, '   +   ', '   -   ')) 
                # Labels = subData$gene_name
                subData$gene_name <- factor(subData$gene_name, levels=subData$gene_name)
                # subData$tissue <- factor(subData$tissue, levels=subData$tissue)

                # subData <- subData %>% 
                # inner_join(tissue_abbr, by = "tissue")
           
                # draw plot 
                p <- ggplot(subData, aes(x=subData$gene_name, y=subData$tissue, size=abs(subData$zscore)))
                p + 
                geom_point(aes(colour=z_score)) + 
                scale_color_manual(values=c('blue', 'brown')) +
                scale_size_continuous(guide=FALSE, range=c(0,max(abs(subData$zscore))))+
                # ggtitle(paste('locus: ', snpsNames[i], '(chromosome', chrosome[i], ')')) +
                labs(x='Gene', y='Tissue') + 
                ggtitle(snpsNames[i])+
                # scale_x_discrete(breaks = subData$gene_name, labels=Labels) + 
                theme(axis.text.x = element_text(size=12, face='bold', angle = 90, hjust = 1)) +
                theme(axis.text.y = element_text(size=12, face='bold')) +
                theme(plot.title = element_text(size=18, face='bold')) +
                # theme(axis.title= element_text(size=18, face='bold')) +
                theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
                theme(legend.position = "none") +
                theme(plot.title = element_text(hjust = 0.5)) + 
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

                # save plot 
                ggsave(paste(snpsNames[i], '_bubble_plot', '.pdf', sep=''), width=12, height=12)
            } 
        """)
        os.chdir(self.srcPath)

        msg = "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done for Bubble Plots!"
        self.logger.info(msg)
        print(msg)

    #  Region Plot 
    def createRegionPlot(self):
        projectName = self.project 

        # current path 
        self.srcPath = os.getcwd()
        self.currentPath = os.getcwd()[:-((len(os.getcwd().split('/')[-1])) + 1)]

        # create directory that holds top gene list (no associated snp)
        self.region_plot_output_path = self.currentPath + '/output/' + projectName  + '/region_plot/'
        os.system("mkdir " + self.region_plot_output_path)

        os.chdir(self.currentPath + '/output/' + projectName + '/top_genes/')
        self.gwas_lead_snp = r("""
              gwas_lead_snp <- read.csv('sorted_top_genes_all_tissues.csv') %>% 
              select(gene_name, chr, start) %>%
              distinct()
            colnames(gwas_lead_snp) = c('snpsNames', 'chrosome', 'startSites')

            """)
        robjects.globalenv['dataframe'] = self.gwas_lead_snp

        # get merged data 
        os.chdir(self.currentPath + '/output/' + projectName + '/annotated_output_files/')
        self.data = r("data <- read.csv('merged_annotated.csv')")
        robjects.globalenv['dataframe'] = self.data 
        os.chdir(self.region_plot_output_path)

        r("""
            # read dataframe 
            startSites = as.numeric(as.character(gwas_lead_snp$startSites))
            snpsNames = as.character(gwas_lead_snp$snpsNames)
            chrosome = as.numeric(as.character(gwas_lead_snp$chrosome))

            # loop through snps site +/- 1000,000 bps 
            for (i in 1:length(startSites)) 
            {
            print(paste("REGION PLOTING", ": ", snpsNames[i]))

            # subset data 
            subData <- subset(data, data$start > startSites[i] - 1000000  & data$start < startSites[i] 
                + 1000000 & data$chr==chrosome[i])
            subData$logp <- -log10(subData$pvalue)
            subData <- subData[order(subData$start),]

            subData <- mutate(subData, sig=ifelse(subData$logp > -log10(0.05/nrow(data)), 'Most Sig', 
            ifelse(subData$logp > 5.30103 & subData$logp <= -log10(0.05/nrow(data)), 'Sig', 
            ifelse(subData$logp > -log10(0.05/nrow(subData)) 
                & subData$logp <= 5.30103, 'Less Sig','Not Sig')))) 
            subData$gene_name <- factor(subData$gene_name, levels=subData$gene_name)

            x_title <- 'Gene'
            y_title <- expression(bold('-log'[10]*'(pvalue)'))

            p <- ggplot(subData, aes(x=gene_name, y=logp))
            p + geom_point(aes(colour = sig)) + 
            scale_color_manual(guide=FALSE, values=c('black', 'black', 'black', 'black')) +
            labs(x=x_title, y=y_title) +
            ggtitle(snpsNames[i])+
            geom_hline(yintercept = 5.30103, linetype='dashed', color='black') +
            geom_hline(yintercept = -log10(0.05/nrow(data)), linetype='solid', color='black') +
            geom_hline(yintercept = -log10(0.05/nrow(subData)), linetype='dotdash', color='black') +
            theme(axis.text.x = element_text(size=12, face='bold', angle = 90, hjust = 1)) +
            theme(axis.text.y = element_text(size=12, face='bold')) + 
            theme(plot.title = element_text(size=18, face='bold')) +
            theme(axis.title.x = element_blank(), axis.title.y=element_text(size=18, face='bold')) + 
            theme(legend.position = "none")+
            theme(plot.title = element_text(hjust = 0.5)) + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

            # output plot 
            ggsave(paste(snpsNames[i], '_region_plot', '.pdf',sep=''), width=12, height=12)
            } 
        """) 

        os.chdir(self.srcPath)

        msg = "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done for Region Plots!"
        self.logger.info(msg)
        print(msg)

    #  locuszoom Plot 
    def createLocuszoomPlot(self):
        projectName = self.project 

        # create a folder to hold locuszoom results 
        locuszoom_plot_path =self.currentPath + '/locuszoom/locuszoom_plots/'
        os.system('mkdir ' + locuszoom_plot_path)

        # get gene lists 
        os.chdir(self.currentPath + '/output/' + projectName + '/top_genes/')
        self.gwas_lead_snp = r("""
              gwas_lead_snp <- read.csv('sorted_top_genes_all_tissues.csv')  %>% 
              select(gene_name, chr, start, end) %>%
              distinct() %>% 
              mutate(flank='1.0MB', run = 'yes', m2zargs= "showAnnot=F")
              # gwas_lead_snp <- gwas_lead_snp[sample(1:nrow(gwas_lead_snp)),]
              colnames(gwas_lead_snp) = c('snp', 'chr', 'start', 'stop', 'flank', 'run', 'm2zargs')
              write.table(gwas_lead_snp, file="batch_locuszoom.txt", sep=" ", quote=FALSE, row.names=F)

            """)
        robjects.globalenv['dataframe'] = self.gwas_lead_snp
        os.chdir(self.srcPath)


        # source files from input path 
        # copy the files including plink, run_locuszoom.py and two other .txt from input path into locuszoom program 
        destination = self.currentPath + '/locuszoom/locuszoom_plots/'     # locuszoom plot destination path 
        plink_destination = self.currentPath + '/locuszoom'


        source = os.listdir(self.top_genes_output_path)

        for file in source: 
            if file.endswith(".txt"):
                shutil.copy(self.top_genes_output_path + file, destination)
                msg =  "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "COPYING: the file '%s' into the folder '%s'" % (file, destination) 
                self.logger.info(msg)
                print(msg)
        os.system('rm ' + self.top_genes_output_path + 'batch_locuszoom.txt')


        # input path 
        self.input_path = self.tools_folder
        source = os.listdir(self.input_path)
        plink_file = self.input_path + 'plink'
        shutil.copy(plink_file, plink_destination)

        msg = "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "COPYING: the file '%s' into the folder '%s'" % (plink_file, plink_destination)
        self.logger.info(msg)
        print(msg)

        for file in source: 
            file = file.strip()
            if file.endswith(".txt"):
                shutil.copy(self.input_path + '/' + file, destination)
                msg = "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "COPYING: the file '%s' into the folder '%s'" % (file, destination)
                self.logger.info(msg) 
                print(msg)

        locuszoom_scrip_cmd = 'run_locuszoom.py'
        shutil.copy(locuszoom_scrip_cmd, destination)
        msg = "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "COPYING: the file '%s' into the folder '%s'" % (locuszoom_scrip_cmd, destination)  
        self.logger.info(msg)
        print(msg)

        # run locuszoom program 
        os.chdir(locuszoom_plot_path) 
        os.system('./run_locuszoom.py')
        msg = "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + 'RUNNING: run_locuszoom.py...'
        self.logger.info(msg)
        print(msg)
        os.chdir(self.currentPath)

        # setup locuszoom plot output file path and move all files to this new output folder 
        locuszoom_plot_files_path = self.currentPath + '/output/' + projectName  + '/locuszoom_plot/'
        os.system("mkdir " + locuszoom_plot_files_path)
        
        src = destination
        dst = locuszoom_plot_files_path
        os.system('mv %s %s' % (src, dst))
        msg = "\n " + datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "MOVING: all files from the folder ' %s ' into the folder ' %s'" % (src, dst)
        self.logger.info(msg)
        print(msg)

        subdirectories = os.listdir(locuszoom_plot_files_path + 'locuszoom_plots/')

        for child in subdirectories: 
            child_path = locuszoom_plot_files_path + 'locuszoom_plots/' + child
            if os.path.isdir(child_path):
               for locuszoom in os.listdir(child_path):
                  if locuszoom.endswith(".pdf"):
                     shutil.copy(child_path + '/' + locuszoom, locuszoom_plot_files_path)
                     os.rename(child_path + '/' + locuszoom, locuszoom_plot_files_path + child.split("_")[-1] + '.pdf')
                     os.remove(locuszoom_plot_files_path + locuszoom)
        path = locuszoom_plot_files_path + 'locuszoom_plots/'
        os.system("rm -rf %s" % path)

        msg = "\n " +  datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "MOVING: all files from the folder '%s' into the folder '%s'" % (src, locuszoom_plot_files_path)
        self.logger.info(msg)
        print(msg)

        msg =  "\n " + datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done for Locuszoom Plots!"
        self.logger.info(msg)
        print(msg)

        msg = "\nElapsed Time: " + timeString(time.time() - self.start_time) + "\n" # calculate how long the program is running
        self.logger.info(msg)
        print(msg)        
        
# --------------------------------------
# main functions 
# --------------------------------------
def main():
    # Instantiate, flush and setup pipeline 
    metaXcanPostprocess = MetaXcanPostprocess()
    metaXcanPostprocess.flushPipeline() 
    metaXcanPostprocess.get_args()
    metaXcanPostprocess.get_parameters()

    # Part One: Annotation 
    metaXcanPostprocess.annotateMetaxcanOutput()

    # Part Two: Top Gene List Without SNPs 
    metaXcanPostprocess.getTopGeneList()
    
    # Part Three: Top Gene List With SNPs  
    if metaXcanPostprocess.multiple_tissue == 'false': 
        metaXcanPostprocess.getTopGeneListWithSNPs()

    # Part Four: QQ-Plot 
    metaXcanPostprocess.createQQPlot()

    # Part Five: Manhattan Plot 
    metaXcanPostprocess.createManhattanPlot()

    # Part Six:  Bubble Plots   
    if metaXcanPostprocess.multiple_tissue == 'false': 
        metaXcanPostprocess.createBubblePlot()
    
    # Part Seven: Region Plots  
    metaXcanPostprocess.createRegionPlot()

	# Part Eight: Locuszoom Plots 
    if metaXcanPostprocess.locuszoom == 'true': 
	    metaXcanPostprocess.createLocuszoomPlot()


# initialize the script
if __name__ == '__main__':
    sys.exit(main())

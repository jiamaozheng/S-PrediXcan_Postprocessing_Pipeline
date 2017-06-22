#!/usr/bin/env python

from helpers import * 

__author__ = "Jiamao Zheng <jiamaoz@yahoo.com>"
__version__ = "Revision: 0.01"
__date__ = "Date: 2017-06-28"

class MetaXcanPostprocess(object):

    def __init__(self, projectName):
        # source code file path 
        self.currentPath = ''       # /src/ path 

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

        # logger 
        self.logger = logging.getLogger()

        # time 
        self.start_time = time.time() 

        # prefix for querying db 
        self.SQL_QUERY_PREFIX = "select e.genename, w.rsid from weights w join extra e on w.gene = e.gene where e.genename = '"
        self.SQL_QUERY_PREFIX_DNG = "select gene, rsid from weights where gene = '"

        # project INFO 
        self.project_time = datetime.now().strftime('%Y-%m-%d') 
        self.project_id = str(myuuid.uuid4())
        # self.project_id = '5e556ffe-2315-40e8-9d88-11f2e61715cb'
        self.project = projectName + '_' + self.project_id + '_' + self.project_time


    # Annotations 
    def annotateMetaxcanOutput(self, projectName):
        projectName = self.project 
        # set up logger  
        log_path = '../log/' + projectName + '/'
        os.makedirs(log_path)
        getLog(self.logger, log_path + projectName + '.log')

        # create directory that holds output and inputs 
        self.output_annotation_path = '../output/' + projectName  + '/annotated_output_files/' 
        os.makedirs(self.output_annotation_path)

        # self.input_path = '../input/'
        self.input_path = '../input/' + projectName +'/'
        os.makedirs(self.input_path)
        os.system('mv ../input/*.* ' + self.input_path)
        os.system('mv ../input/plink ' + self.input_path)


        # get input file lists (raw metaxcan output files *.csv) 
        inputFileList = glob.glob(self.input_path + '/*.csv')

        # loop through file lists 
        for inputFilename in inputFileList:
            msg = "ANNOTATING GENES: " + inputFilename 
            self.logger.info(msg)
            print(msg)
            # read data 
            r("data <- read.csv('%s')" %inputFilename)
            data = r('na.omit(data)')
            robjects.globalenv['dataframe'] = data

            # annotation library 
            grch37 = r('grch37')
            robjects.globalenv['dataframe'] = grch37

            # annotating  
            # if inputFilename == 'xxxxxDGN-WB-unscaled.csv':
            #     annotatedData = r("inner_join(data[, c('gene', 'gene_name', 'zscore', 'effect_size', 'pvalue', 'n_snps_in_model', 'pred_perf_p', 'pred_perf_R2')], grch37[, c('symbol', 'chr', 'start', 'end')], by=c('gene'='symbol'))")
            # else:
            annotatedData = r("inner_join(data[, c('gene', 'gene_name', 'zscore', 'effect_size', 'pvalue', 'n_snps_in_model', 'pred_perf_pval', 'pred_perf_r2')], grch37[, c('ensgene', 'chr', 'start', 'end')], by =c('gene'='ensgene'))")
            annotatedData = annotatedData.drop_duplicates()

            # ouput annotated data 
            annotatedData.to_csv(self.output_annotation_path + inputFilename.split('/')[-1][:-4] + '_annotated.csv', index=None)

        # merge all annotated metaxcan files 
        # get annotated file lists 
        self.outputFileListWithoutSNPs = glob.glob(self.output_annotation_path + "*.csv")

        dfListWithoutSNPs = []
        for outputFileName in self.outputFileListWithoutSNPs:
            msg = "CONCATING FILE: " + outputFileName 
            self.logger.info(msg)
            print(msg)

            r("data <- read.csv('%s')" %outputFileName)
            data = r('data <- na.omit(data)')
            robjects.globalenv['dataframe'] = data

            if 'CrossTissue_elasticNet' in outputFileName:
                outputFileName = outputFileName[:-25]
            elif 'DGN-WB-unscaled' in outputFileName: 
                outputFileName = outputFileName [:-23]
            else: 
                outputFileName = outputFileName[:-25]
                outputFileName = outputFileName[3:]

            tissue = outputFileName.split('/')[-1][3:]
            data.insert(2, 'tissue', tissue)
            dfListWithoutSNPs.append(data)
            
        concatDfWithoutSNPs = pandas.concat(dfListWithoutSNPs, axis = 0)

        #output merged data 
        self.merged_output = self.output_annotation_path + 'merged_output.csv'
        concatDfWithoutSNPs.to_csv(self.merged_output, index = None)
        msg = datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done for Annotating metaxcan output!"
        self.logger.info(msg)
        print(msg)

    #  QQ-Plot  
    def createQQPlot(self, projectName):
        projectName = self.project 

        # create directory that holds qq plot output
        self.output_qqplot_path = '../output/' + projectName  + '/qq-plot/' 
        os.system("mkdir " + self.output_qqplot_path)

        # loop through file lists 
        for outputFileName in self.outputFileListWithoutSNPs: 
            msg = "QQ-PLOT(no snps): " + outputFileName 
            self.logger.info(msg)
            print(msg)
            # get GWAS output data 
            r("data <- read.csv('%s')" %outputFileName)
            data = r('data <- na.omit(data)')
            robjects.globalenv['dataframe'] = data

            # draw qq-plot and save them to files 
            r.pdf('%s%s%s%s'%(self.output_qqplot_path, 'qq-plot_', outputFileName.split('/')[-1][:-4],'.pdf'))
            qqman.qq(data['pvalue'],  main = "")
            r['dev.off']()
        msg = "QQ-PLOT(no snps): " + "merged_output"
        self.logger.info(msg)
        print(msg)

        # set up path and read merged annotated metaxcan output files  
        r("data <- read.csv('%s')" % (self.merged_output))
        data = r('data <- na.omit(data)')
        robjects.globalenv['dataframe'] = data

        # draw qq-plot and save them to files 
        r.pdf('%s%s%s%s'%(self.output_qqplot_path,'qq-plot_', 'merged_output', '.pdf'))
        qqman.qq(data['pvalue'],  main = "QQ plot")
        r['dev.off']()

        msg = datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done for QQ Plot!"
        self.logger.info(msg)
        print(msg)
    
    #  Manhattan-Plot 
    def createManhattanPlot(self, projectName):
        projectName = self.project 

        # create directory that holds manhattan plot output
        self.manhattanplot_output_path = '../output/' + projectName  + '/manhattan-plot/'
        os.system("mkdir " + self.manhattanplot_output_path)

        # loop through file lists 
        for outputFileName in self.outputFileListWithoutSNPs: 
            msg = "MANHATTAN-PLOT(no snps): " + outputFileName
            self.logger.info(msg)
            print(msg)
            # read data 
            r("data <- read.csv('%s')" %outputFileName)
            r('data$chr <- as.numeric(as.character(data$chr))')
            data = r('data <- na.omit(data)')
            robjects.globalenv['dataframe'] = data

            # draw manhattan and save them to files 
            r.pdf('%s%s%s%s'%(self.manhattanplot_output_path, 'manhattan-Plot_', outputFileName.split('/')[-1][:-4], '.pdf'))
            qqman.manhattan(data, chr = 'chr', bp='start', p='pvalue', snp='gene_name', main = 'Manhattan plot')
            r['dev.off']()

        # read data 
        r("data <- read.csv('%s')" % self.merged_output)
        r('data$chr <- as.numeric(as.character(data$chr))')
        data = r('data <- na.omit(data)')
        robjects.globalenv['dataframe'] = data

        # draw manhattan and save them to files 
        r.pdf('%s%s%s%s'%(self.manhattanplot_output_path, 'manhattan-plot_', 'merged_output', '.pdf'))
        qqman.manhattan(data, chr = 'chr', bp='start', p='pvalue', snp='gene_name', main = 'Manhattan plot')
        r['dev.off']()

        msg = datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done for Manhattan Plot!"
        self.logger.info(msg)
        print(msg)

    # Top Gene List without SNPs
    def getTopGeneList(self, projectName):
        projectName = self.project 

        # create directory that holds top gene list (no associated snp)
        self.top_genes_output_path = '../output/' + projectName  +'/top_genes/'
        os.system("mkdir " + self.top_genes_output_path)

        # loop through files 
        for filename in self.outputFileListWithoutSNPs:
            msg = "TOP GENE LIST: " + filename
            self.logger.info(msg)
            print(msg)

            df = pandas.read_csv(filename) 
            if 'CrossTissue_elasticNet' in filename:
                tissue_name = filename[:-25]
            elif 'DGN-WB-unscaled' in filename: 
                tissue_name = filename [:-23]
            else: 
                tissue_name = filename[:-25]
                tissue_name = tissue_name[3:]
            df.insert(2, 'tissue', tissue_name.split('/')[-1])

            # sort data by defined column 
            df.sort_values(['pvalue', 'n_snps_in_model'], ascending=[1, 0], inplace=True)
            total_rows = df.shape[0] # number of row count shape[1] is number of col count 

            # total_rows_in_total += total_rows 
            cut_off = 0.05/total_rows 

            # cut_off_list.append(cut_off)
            top_gene_list = df[df['pvalue'] < cut_off]
            # top_gene_list = top_gene_list[top_gene_list['pred_perf_R2'] > 0.01]
            # top_gene_list.drop('gene', axis=1, inplace=True)
            msg = 'CALCULATING TISSUE-WIDE P-VALUE: ' + str(cut_off)
            self.logger.info(msg)
            print(msg)

            #output data 
            top_gene_list.to_csv(self.top_genes_output_path + "sorted_top_genes_%s"%tissue_name.split('/')[-1] + ".csv", index = None)

            # dfList.append(top_gene_list)

        df = pandas.read_csv(self.merged_output) 
        total_rows_in_total = df.shape[0] 
        msg = 'total_row_in_total: ' + str(total_rows_in_total)
        self.logger.info(msg)
        print(msg)

        total_cut_off = 0.05/total_rows_in_total 
        msg = 'CALCULATING GENOME-WIDE P-VALUE: ' + str(total_cut_off)
        self.logger.info(msg)
        print(msg)

        df = df[df['pvalue'] < total_cut_off]

        # sort data by defined column 
        df.sort_values(['pvalue', 'n_snps_in_model'], ascending=[1, 0], inplace=True)

        #output data 
        self.sorted_top_genes_all_tissues_path = self.top_genes_output_path+ 'sorted_top_genes_all_tissues.csv' 
        df.to_csv(self.sorted_top_genes_all_tissues_path, index = None)
        self.top_gene_list_dataframe = df

        msg = datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done for Top Gene List!"
        self.logger.info(msg)
        print(msg)

    #  Top Gene List with SNPs 
    def getTopGeneListWithSNPs(self, projectName):
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
        dbFileList = glob.glob(self.input_path + "/*.db")

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

                    msg = 'FETCH SNPs for GENE %s ' % gene_lists[k] +  'FROM DATABASE: %s' % database_names[i]
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

        msg = datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done for Top Gene List with SNPs!"
        self.logger.info(msg)
        print(msg)

    #  Bubble Plot 
    def createBubblePlot(self, projectName):
        projectName = self.project 

        # create directory that holds top gene list (no associated snp)
        self.bubble_plot_output_path = '../output/' + projectName  + '/bubble_plot/'
        os.system("mkdir " + self.bubble_plot_output_path)

        # current path 
        self.currentPath = os.getcwd()
        os.chdir(self.currentPath[:-4] + '/input/' + projectName +'/')

        # get tissue abbr name 
        self.tissue_abbr = r("""

        tissue_abbr <- read.delim('gtex_tissue_abbr.txt', sep='\t') %>% 
        select(tissue_site_detail_abbr,tissue_site_detail_id)
        colnames(tissue_abbr)[colnames(tissue_abbr) == 'tissue_site_detail_id'] = "tissue"

        """) 
        robjects.globalenv['dataframe'] = self.tissue_abbr

        # get gene lists 
        os.chdir(self.currentPath[:-4] + '/output/' + projectName + '/top_genes/')
        self.gwas_lead_snp = r("""
              gwas_lead_snp <- read.csv('sorted_top_genes_all_tissues.csv') %>% 
              select(gene_name, chr, start) %>%
              distinct()
            colnames(gwas_lead_snp) = c('snpsNames', 'chrosome', 'startSites')

            """)
        robjects.globalenv['dataframe'] = self.gwas_lead_snp

        # get merged data 
        os.chdir(self.currentPath[:-4] + '/output/' + projectName + '/annotated_output_files/')
        self.data = r("data <- read.csv('merged_output.csv')")
        robjects.globalenv['dataframe'] = self.data 

        os.chdir(self.currentPath[:-4] + '/output/' + projectName + '/bubble_plot/')

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

                subData <- subData %>% 
                inner_join(tissue_abbr, by = "tissue")
           
                # draw plot 
                p <- ggplot(subData, aes(x=subData$gene_name, y=subData$tissue_site_detail_abbr, size=abs(subData$zscore)))
                p + 
                geom_point(aes(colour=z_score)) + 
                scale_color_manual(values=c('blue', 'brown')) +
                scale_size_continuous(guide=FALSE, range=c(0,max(abs(subData$zscore))))+
                # ggtitle(paste('locus: ', snpsNames[i], '(chromosome', chrosome[i], ')')) +
                labs(x='Gene', y='Tissue') + 
                ggtitle(snpsNames[i])+
                # scale_x_discrete(breaks = subData$gene_name, labels=Labels) + 
                theme(axis.text.x = element_text(size=16, face='bold', angle = 90, hjust = 1)) +
                theme(axis.text.y = element_text(size=16, face='bold')) +
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
        os.chdir(self.currentPath)

        msg = datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done for Bubble Plots!"
        self.logger.info(msg)
        print(msg)

    #  Region Plot 
    def createRegionPlot(self, projectName):
        projectName = self.project 

        # create directory that holds top gene list (no associated snp)
        self.region_plot_output_path = '../output/' + projectName  + '/region_plot/'
        os.system("mkdir " + self.region_plot_output_path)

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

            p <- ggplot(subData, aes(x=gene_name, y=logp))
            p + geom_point(aes(colour = sig)) + 
            scale_color_manual(guide=FALSE, values=c('black', 'black', 'black', 'black')) +
                # ggtitle(paste("locus: ", snpsNames[i], '(chromosome', chrosome[i], ')')) + 
            labs(x='Gene', y='-log10(p-value)') +
            ggtitle(snpsNames[i])+
            geom_hline(yintercept = 5.30103, linetype='dashed', color='black') +
            geom_hline(yintercept = -log10(0.05/nrow(data)), color='black') +
            geom_hline(yintercept = -log10(0.05/nrow(subData)), linetype='dashed', color='black') +
            theme(axis.text.x = element_text(size=16, face='bold', angle = 90, hjust = 1)) +
            theme(axis.text.y = element_text(size=16, face='bold')) + 
            theme(plot.title = element_text(size=18, face='bold')) +
            theme(axis.title.x = element_blank(), axis.title.y=element_text(size=18, face='bold')) + 
            theme(legend.position = "none")+
            theme(plot.title = element_text(hjust = 0.5)) + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

            # output plot 
            ggsave(paste(snpsNames[i], '_region_plot', '.pdf',sep=''), width=12, height=12)
            } 
        """) 

        os.chdir(self.currentPath)

        msg = datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done for Region Plots!"
        self.logger.info(msg)
        print(msg)

    #  locuszoom Plot 
    def createLocuszoomPlot(self, projectName):
        projectName = self.project 

        # create a folder to hold locuszoom results 
        locuszoom_plot_path = '../locuszoom/locuszoom_plots/'
        os.system('mkdir ' + locuszoom_plot_path)

        # get gene lists 
        os.chdir(self.currentPath[:-4] + '/output/' + projectName + '/top_genes/')
        self.gwas_lead_snp = r("""
              gwas_lead_snp <- read.csv('sorted_top_genes_all_tissues.csv')  %>% 
              select(gene_name, chr, start, end) %>%
              distinct() %>% 
              mutate(flank='1.0MB', run = 'yes', m2zargs= "showAnnot=F")
              colnames(gwas_lead_snp) = c('snp', 'chr', 'start', 'stop', 'flank', 'run', 'm2zargs')
              write.table(gwas_lead_snp, file="batch_locuszoom.txt", sep=" ", quote=FALSE, row.names=F)

            """)
        robjects.globalenv['dataframe'] = self.gwas_lead_snp
        os.chdir(self.currentPath)


        # source files from input path 
        # copy the files including plink, run_locuszoom.py and two other .txt from input path into locuszoom program 
        destination = self.currentPath[:-4] + '/locuszoom/locuszoom_plots/'     # locuszoom plot destination path 
        plink_destination = self.currentPath[:-4] + '/locuszoom'


        source = os.listdir(self.top_genes_output_path)

        for file in source: 
            if file.endswith(".txt"):
                shutil.copy(self.top_genes_output_path + file, destination)
                msg = "COPYING: the file '%s' into the folder '%s'" % (file, destination) 
                self.logger.info(msg)
                print(msg)
        os.system('rm ' + self.top_genes_output_path + 'batch_locuszoom.txt')


        # input path 
        source = os.listdir(self.input_path)

        plink_file = self.input_path + '/plink'
        shutil.copy(plink_file, plink_destination)

        msg = "COPYING: the file '%s' into the folder '%s'" % (plink_file, plink_destination)
        self.logger.info(msg)
        print(msg)

        for file in source: 
            file = file.strip()
            if file.endswith(".txt"):
                shutil.copy(self.input_path + '/' + file, destination)
                msg = "COPYING: the file '%s' into the folder '%s'" % (file, destination)
                self.logger.info(msg) 
                print(msg)

        locuszoom_scrip_cmd = 'run_locuszoom.py'
        shutil.copy(locuszoom_scrip_cmd, destination)
        msg = "COPYING: the file '%s' into the folder '%s'" % (locuszoom_scrip_cmd, destination)  
        self.logger.info(msg)
        print(msg)

        # run locuszoom program 
        os.chdir(locuszoom_plot_path) 
        os.system('./run_locuszoom.py')
        msg = 'RUNNING: run_locuszoom.py...'
        self.logger.info(msg)
        print(msg)
        os.chdir(self.currentPath)

        # setup locuszoom plot output file path and move all files to this new output folder 
        locuszoom_plot_files_path = '../output/' + projectName  + '/locuszoom_plot/'
        os.system("mkdir " + locuszoom_plot_files_path)
        
        src = destination
        dst = locuszoom_plot_files_path
        os.system('mv %s %s' % (src, dst))
        msg = "MOVING: all files from the folder ' %s ' into the folder ' %s'" % (src, dst)
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

        msg = "MOVING: all files from the folder '%s' into the folder '%s'" % (src, locuszoom_plot_files_path)
        self.logger.info(msg)
        print(msg)

        msg = datetime.now().strftime('%Y.%m.%d.%H:%M:%S ') + "Done for Locuszoom Plots!"
        self.logger.info(msg)
        print(msg)

        msg = "\nElapsed Time: " + timeString(time.time() - self.start_time) + "\n" # calculate how long the program is running
        self.logger.info(msg)
        print(msg)
        
# --------------------------------------
# main functions 
# --------------------------------------
def main():
    # Instantiate Class ------------------------------------
    metaXcanPostprocess = MetaXcanPostprocess(sys.argv[1]) 

    # Part One: Annotation ------------------------------------
    metaXcanPostprocess.annotateMetaxcanOutput(sys.argv[1])

    # Part Two: QQ-Plot ------------------------------------
    metaXcanPostprocess.createQQPlot(sys.argv[1]) 

    # Part Three: Manhattan Plot ------------------------------------
    metaXcanPostprocess.createManhattanPlot(sys.argv[1])

    # Part Four: Top Gene List Without SNPs ------------------------------------
    metaXcanPostprocess.getTopGeneList(sys.argv[1])

    # Part Five: Top Gene List With SNPs  ------------------------------------
    metaXcanPostprocess.getTopGeneListWithSNPs(sys.argv[1])

    #  Part Six:  Bubble Plots   ------------------------------------
    metaXcanPostprocess.createBubblePlot(sys.argv[1])
    
    # Part Seven: Region Plots  ------------------------------------
    metaXcanPostprocess.createRegionPlot(sys.argv[1])

    # Part Eight: Locuszoom Plots ------------------------------------
    metaXcanPostprocess.createLocuszoomPlot(sys.argv[1])



# initialize the script
if __name__ == '__main__':
    sys.exit(main())





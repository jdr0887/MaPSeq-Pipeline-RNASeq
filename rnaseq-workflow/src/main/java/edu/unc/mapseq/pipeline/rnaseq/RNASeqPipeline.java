package edu.unc.mapseq.pipeline.rnaseq;

import java.io.File;
import java.util.HashSet;
import java.util.List;
import java.util.ResourceBundle;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.renci.jlrm.condor.CondorJob;
import org.renci.jlrm.condor.CondorJobEdge;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.dao.MaPSeqDAOException;
import edu.unc.mapseq.dao.model.HTSFSample;
import edu.unc.mapseq.dao.model.SequencerRun;
import edu.unc.mapseq.module.bedtools.CoverageBedCLI;
import edu.unc.mapseq.module.core.GUnZipCLI;
import edu.unc.mapseq.module.core.MoveCLI;
import edu.unc.mapseq.module.core.RemoveCLI;
import edu.unc.mapseq.module.filter.PruneISOFormsFromGeneQuantFileCLI;
import edu.unc.mapseq.module.filter.StripTrailingTabsCLI;
import edu.unc.mapseq.module.mapsplice.MapSpliceCLI;
import edu.unc.mapseq.module.mapsplice.RSEMCalculateExpressionCLI;
import edu.unc.mapseq.module.picard.PicardAddOrReplaceReadGroupsCLI;
import edu.unc.mapseq.module.picard.PicardSortOrderType;
import edu.unc.mapseq.module.qc.NormBedExonQuantCLI;
import edu.unc.mapseq.module.qc.NormalizeQuartileCLI;
import edu.unc.mapseq.module.samtools.SAMToolsFlagstatCLI;
import edu.unc.mapseq.module.samtools.SAMToolsIndexCLI;
import edu.unc.mapseq.module.samtools.SAMToolsSortCLI;
import edu.unc.mapseq.module.samtools.SortByReferenceAndNameCLI;
import edu.unc.mapseq.module.ubu.UBUFastqFormatterCLI;
import edu.unc.mapseq.module.ubu.UBUSamFilterCLI;
import edu.unc.mapseq.module.ubu.UBUSamJunctionCLI;
import edu.unc.mapseq.module.ubu.UBUSamTranslateCLI;
import edu.unc.mapseq.pipeline.AbstractPipeline;
import edu.unc.mapseq.pipeline.PipelineException;
import edu.unc.mapseq.pipeline.PipelineJobFactory;
import edu.unc.mapseq.pipeline.PipelineUtil;

public class RNASeqPipeline extends AbstractPipeline<RNASeqPipelineBeanService> {

    private final Logger logger = LoggerFactory.getLogger(RNASeqPipeline.class);

    private RNASeqPipelineBeanService pipelineBeanService;

    public RNASeqPipeline() {
        super();
    }

    @Override
    public String getName() {
        return RNASeqPipeline.class.getSimpleName().replace("Pipeline", "");
    }

    @Override
    public String getVersion() {
        ResourceBundle bundle = ResourceBundle.getBundle("edu/unc/mapseq/pipeline/rnaseq/pipeline");
        String version = bundle.getString("version");
        return StringUtils.isNotEmpty(version) ? version : "0.0.1-SNAPSHOT";
    }

    @Override
    public Graph<CondorJob, CondorJobEdge> createGraph() throws PipelineException {
        logger.info("ENTERING createGraph()");

        DirectedGraph<CondorJob, CondorJobEdge> graph = new DefaultDirectedGraph<CondorJob, CondorJobEdge>(
                CondorJobEdge.class);

        int count = 0;

        if (getWorkflowPlan().getSequencerRun() == null && getWorkflowPlan().getHTSFSamples() == null) {
            logger.error("Don't have either sequencerRun and htsfSample");
            throw new PipelineException("Don't have either sequencerRun and htsfSample");
        }

        Set<HTSFSample> htsfSampleSet = new HashSet<HTSFSample>();

        if (getWorkflowPlan().getSequencerRun() != null) {
            logger.info("sequencerRun: {}", getWorkflowPlan().getSequencerRun().toString());
            try {
                htsfSampleSet.addAll(this.pipelineBeanService.getMaPSeqDAOBean().getHTSFSampleDAO()
                        .findBySequencerRunId(getWorkflowPlan().getSequencerRun().getId()));
            } catch (MaPSeqDAOException e) {
                e.printStackTrace();
            }
        }

        if (getWorkflowPlan().getHTSFSamples() != null) {
            htsfSampleSet.addAll(getWorkflowPlan().getHTSFSamples());
        }

        logger.info("htsfSampleSet.size(): {}", htsfSampleSet.size());

        for (HTSFSample htsfSample : htsfSampleSet) {

            if ("Undetermined".equals(htsfSample.getBarcode())) {
                continue;
            }

            SequencerRun sequencerRun = htsfSample.getSequencerRun();
            File outputDirectory = createOutputDirectory(sequencerRun.getName(), htsfSample.getName(), getName());
            File tmpDir = new File(outputDirectory, "tmp");
            tmpDir.mkdirs();

            List<File> readPairList = PipelineUtil.getReadPairList(htsfSample.getFileDatas(), sequencerRun.getName(),
                    htsfSample.getLaneIndex());
            logger.info("fileList = {}", readPairList.size());

            // assumption: a dash is used as a delimiter between a participantId
            // and the external code
            int idx = htsfSample.getName().lastIndexOf("-");
            String participantId = idx != -1 ? htsfSample.getName().substring(0, idx) : htsfSample.getName();

            if (readPairList.size() == 2) {

                File r1FastqFile = readPairList.get(0);
                String r1FastqRootName = PipelineUtil.getRootFastqName(r1FastqFile.getName());

                File r2FastqFile = readPairList.get(1);
                String r2FastqRootName = PipelineUtil.getRootFastqName(r2FastqFile.getName());

                String fastqLaneRootName = StringUtils.removeEnd(r2FastqRootName, "_R2");

                try {
                    // new job
                    CondorJob gunzipFastqR1Job = PipelineJobFactory.createJob(++count, GUnZipCLI.class,
                            getWorkflowPlan(), htsfSample);
                    gunzipFastqR1Job.addArgument(GUnZipCLI.GZFILE, r1FastqFile.getAbsolutePath());
                    File gunzippedFastqR1 = new File(outputDirectory, r1FastqRootName + ".fastq");
                    gunzipFastqR1Job.addArgument(GUnZipCLI.EXTRACTFILE, gunzippedFastqR1.getAbsolutePath());
                    graph.addVertex(gunzipFastqR1Job);

                    // new job
                    CondorJob fastqFormatterR1Job = PipelineJobFactory.createJob(++count, UBUFastqFormatterCLI.class,
                            getWorkflowPlan(), htsfSample);
                    fastqFormatterR1Job.addArgument(UBUFastqFormatterCLI.INPUT, gunzippedFastqR1.getAbsolutePath());
                    fastqFormatterR1Job.addArgument(UBUFastqFormatterCLI.STRIP);
                    fastqFormatterR1Job.addArgument(UBUFastqFormatterCLI.SUFFIX, "/1");
                    File fastqFormatterR1Out = new File(outputDirectory, gunzippedFastqR1.getName().replace(".fastq",
                            ".filtered.fastq"));
                    fastqFormatterR1Job.addArgument(UBUFastqFormatterCLI.OUTPUT, fastqFormatterR1Out.getAbsolutePath());
                    graph.addVertex(fastqFormatterR1Job);
                    graph.addEdge(gunzipFastqR1Job, fastqFormatterR1Job);

                    // new job
                    CondorJob removeR1Job = PipelineJobFactory.createJob(++count, RemoveCLI.class, getWorkflowPlan(),
                            htsfSample);
                    removeR1Job.addArgument(RemoveCLI.FILE, gunzippedFastqR1.getAbsolutePath());
                    graph.addVertex(removeR1Job);
                    graph.addEdge(fastqFormatterR1Job, removeR1Job);

                    // new job
                    CondorJob gunzipFastqR2Job = PipelineJobFactory.createJob(++count, GUnZipCLI.class,
                            getWorkflowPlan(), htsfSample);
                    gunzipFastqR2Job.addArgument(GUnZipCLI.GZFILE, r2FastqFile.getAbsolutePath());
                    File gunzippedFastqR2 = new File(outputDirectory, r2FastqRootName + ".fastq");
                    gunzipFastqR2Job.addArgument(GUnZipCLI.EXTRACTFILE, gunzippedFastqR2.getAbsolutePath());
                    graph.addVertex(gunzipFastqR2Job);

                    // new job
                    CondorJob fastqFormatterR2Job = PipelineJobFactory.createJob(++count, UBUFastqFormatterCLI.class,
                            getWorkflowPlan(), htsfSample);
                    fastqFormatterR2Job.addArgument(UBUFastqFormatterCLI.INPUT, gunzippedFastqR2.getAbsolutePath());
                    fastqFormatterR2Job.addArgument(UBUFastqFormatterCLI.STRIP);
                    fastqFormatterR2Job.addArgument(UBUFastqFormatterCLI.SUFFIX, "/2");
                    File fastqFormatterR2Out = new File(outputDirectory, gunzippedFastqR2.getName().replace(".fastq",
                            ".filtered.fastq"));
                    fastqFormatterR2Job.addArgument(UBUFastqFormatterCLI.OUTPUT, fastqFormatterR2Out.getAbsolutePath());
                    graph.addVertex(fastqFormatterR2Job);
                    graph.addEdge(gunzipFastqR2Job, fastqFormatterR2Job);

                    // new job
                    CondorJob removeR2Job = PipelineJobFactory.createJob(++count, RemoveCLI.class, getWorkflowPlan(),
                            htsfSample);
                    removeR2Job.addArgument(RemoveCLI.FILE, gunzippedFastqR2.getAbsolutePath());
                    graph.addVertex(removeR2Job);
                    graph.addEdge(fastqFormatterR2Job, removeR2Job);

                    // new job
                    CondorJob mapspliceJob = PipelineJobFactory.createJob(++count, MapSpliceCLI.class,
                            getWorkflowPlan(), htsfSample);
                    mapspliceJob.addArgument(MapSpliceCLI.BAM);
                    mapspliceJob.addArgument(MapSpliceCLI.FUSIONNONCANONICAL);
                    mapspliceJob.addArgument(MapSpliceCLI.QUALSCALE, "phred33");
                    mapspliceJob
                            .addArgument(MapSpliceCLI.ALLCHROMOSOMEFILES,
                                    "/proj/seq/LBG/tier1data/nextgenseq/seqware-analysis/mapsplice_rsem/Genomes/hg19_M_rCRS/hg19_M_rCRS.fa");
                    mapspliceJob.addArgument(MapSpliceCLI.THREADS, 8);
                    mapspliceJob.setNumberOfProcessors(8);
                    mapspliceJob.addArgument(MapSpliceCLI.BOWTIEINDEXPATH,
                            this.pipelineBeanService.getBowtieIndexDirectory());
                    mapspliceJob.addArgument(MapSpliceCLI.CHROMOSOMEDIRECTORY,
                            this.pipelineBeanService.getChromosomeDirectory());
                    mapspliceJob.addArgument(MapSpliceCLI.FASTQR1, fastqFormatterR1Out.getAbsolutePath());
                    mapspliceJob.addArgument(MapSpliceCLI.FASTQR2, fastqFormatterR2Out.getAbsolutePath());
                    mapspliceJob.addArgument(MapSpliceCLI.OUTPUT, outputDirectory.getAbsolutePath());
                    graph.addVertex(mapspliceJob);
                    graph.addEdge(fastqFormatterR1Job, mapspliceJob);
                    graph.addEdge(fastqFormatterR2Job, mapspliceJob);

                    // new job
                    CondorJob moveJob = PipelineJobFactory.createJob(++count, MoveCLI.class, getWorkflowPlan(),
                            htsfSample, false);
                    File mapspliceOutput = new File(outputDirectory, "alignments.bam");
                    moveJob.addArgument(MoveCLI.SOURCE, mapspliceOutput.getAbsolutePath());
                    File mapspliceRenamedOutput = new File(outputDirectory, fastqLaneRootName + ".bam");
                    moveJob.addArgument(MoveCLI.DESTINATION, mapspliceRenamedOutput.getAbsolutePath());
                    graph.addVertex(moveJob);
                    graph.addEdge(mapspliceJob, moveJob);

                    // new job
                    CondorJob picardAddOrReplaceReadGroupsJob = PipelineJobFactory.createJob(++count,
                            PicardAddOrReplaceReadGroupsCLI.class, getWorkflowPlan(), htsfSample);
                    picardAddOrReplaceReadGroupsJob.addArgument(PicardAddOrReplaceReadGroupsCLI.INPUT,
                            mapspliceRenamedOutput.getAbsolutePath());
                    File fixRGOutput = new File(outputDirectory, mapspliceRenamedOutput.getName().replace(".bam",
                            ".fixed-rg.bam"));
                    picardAddOrReplaceReadGroupsJob.addArgument(PicardAddOrReplaceReadGroupsCLI.OUTPUT,
                            fixRGOutput.getAbsolutePath());
                    picardAddOrReplaceReadGroupsJob.addArgument(PicardAddOrReplaceReadGroupsCLI.TMPDIR,
                            tmpDir.getAbsolutePath());
                    picardAddOrReplaceReadGroupsJob.addArgument(PicardAddOrReplaceReadGroupsCLI.SORTORDER,
                            PicardSortOrderType.COORDINATE.toString().toLowerCase());
                    picardAddOrReplaceReadGroupsJob.addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPID, htsfSample
                            .getId().toString());
                    picardAddOrReplaceReadGroupsJob.addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPLIBRARY,
                            participantId);
                    picardAddOrReplaceReadGroupsJob.addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPPLATFORM,
                            sequencerRun.getPlatform().getInstrument());
                    picardAddOrReplaceReadGroupsJob.addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPPLATFORMUNIT,
                            sequencerRun.getPlatform().getInstrumentModel());
                    picardAddOrReplaceReadGroupsJob.addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPSAMPLENAME,
                            participantId);
                    graph.addVertex(picardAddOrReplaceReadGroupsJob);
                    graph.addEdge(moveJob, picardAddOrReplaceReadGroupsJob);

                    // new job
                    CondorJob samtoolsSortJob = PipelineJobFactory.createJob(++count, SAMToolsSortCLI.class,
                            getWorkflowPlan(), htsfSample);
                    samtoolsSortJob.addArgument(SAMToolsSortCLI.INPUT, fixRGOutput.getAbsolutePath());
                    File samtoolsSortOut = new File(outputDirectory, fixRGOutput.getName().replace(".bam",
                            ".sorted.bam"));
                    samtoolsSortJob.addArgument(SAMToolsSortCLI.OUTPUT, samtoolsSortOut.getAbsolutePath());
                    graph.addVertex(samtoolsSortJob);
                    graph.addEdge(picardAddOrReplaceReadGroupsJob, samtoolsSortJob);

                    // new job
                    CondorJob samtoolsIndexJob = PipelineJobFactory.createJob(++count, SAMToolsIndexCLI.class,
                            getWorkflowPlan(), htsfSample);
                    samtoolsIndexJob.addArgument(SAMToolsIndexCLI.INPUT, samtoolsSortOut.getAbsolutePath());
                    File samtoolsIndexOutput = new File(outputDirectory, samtoolsSortOut.getName().replace(".bam",
                            ".bai"));
                    samtoolsIndexJob.addArgument(SAMToolsIndexCLI.OUTPUT, samtoolsIndexOutput.getAbsolutePath());
                    graph.addVertex(samtoolsIndexJob);
                    graph.addEdge(samtoolsSortJob, samtoolsIndexJob);

                    // new job
                    CondorJob samtoolsFlagstatJob = PipelineJobFactory.createJob(++count, SAMToolsFlagstatCLI.class,
                            getWorkflowPlan(), htsfSample);
                    samtoolsFlagstatJob.addArgument(SAMToolsFlagstatCLI.INPUT, samtoolsSortOut.getAbsolutePath());
                    File samtoolsFlagstatOut = new File(outputDirectory, samtoolsSortOut.getName().replace(".bam",
                            ".flagstat"));
                    samtoolsFlagstatJob.addArgument(SAMToolsFlagstatCLI.OUTPUT, samtoolsFlagstatOut.getAbsolutePath());
                    graph.addVertex(samtoolsFlagstatJob);
                    graph.addEdge(samtoolsIndexJob, samtoolsFlagstatJob);

                    // new job
                    CondorJob ubuSamJunctionJob = PipelineJobFactory.createJob(++count, UBUSamJunctionCLI.class,
                            getWorkflowPlan(), htsfSample);
                    ubuSamJunctionJob.addArgument(UBUSamJunctionCLI.JUNCTIONS, this.pipelineBeanService.getJunctions());
                    ubuSamJunctionJob.addArgument(UBUSamJunctionCLI.INPUT, samtoolsSortOut.getAbsolutePath());
                    File ubuSamJunctionOut = new File(outputDirectory, samtoolsSortOut.getName().replace(".bam",
                            ".junction_quantification.txt"));
                    ubuSamJunctionJob.addArgument(UBUSamJunctionCLI.OUTPUT, ubuSamJunctionOut.getAbsolutePath());
                    graph.addVertex(ubuSamJunctionJob);
                    graph.addEdge(samtoolsIndexJob, ubuSamJunctionJob);

                    // new job
                    CondorJob coverageBedJob = PipelineJobFactory.createJob(++count, CoverageBedCLI.class,
                            getWorkflowPlan(), htsfSample);
                    coverageBedJob.addArgument(CoverageBedCLI.INPUT, samtoolsSortOut.getAbsolutePath());
                    coverageBedJob.addArgument(CoverageBedCLI.BED, this.getPipelineBeanService().getCompositeExons());
                    coverageBedJob.addArgument(CoverageBedCLI.SPLITBED);
                    File coverageBedOut = new File(outputDirectory, samtoolsSortOut.getName().replace(".bam",
                            ".coverageBedOut.txt"));
                    coverageBedJob.addArgument(CoverageBedCLI.OUTPUT, coverageBedOut.getAbsolutePath());
                    graph.addVertex(coverageBedJob);
                    graph.addEdge(samtoolsSortJob, coverageBedJob);

                    // new job
                    CondorJob normBedExonQuantJob = PipelineJobFactory.createJob(++count, NormBedExonQuantCLI.class,
                            getWorkflowPlan(), htsfSample);
                    normBedExonQuantJob.addArgument(NormBedExonQuantCLI.INFILE, coverageBedOut.getAbsolutePath());
                    normBedExonQuantJob.addArgument(NormBedExonQuantCLI.COMPOSITEBED, this.getPipelineBeanService()
                            .getCompositeExons());
                    File normBedExonQuantOut = new File(outputDirectory, coverageBedOut.getName().replace(
                            ".coverageBedOut.txt", ".normBedExonQuantOut.txt"));
                    normBedExonQuantJob.addArgument(NormBedExonQuantCLI.OUTFILE, normBedExonQuantOut.getAbsolutePath());
                    graph.addVertex(normBedExonQuantJob);
                    graph.addEdge(coverageBedJob, normBedExonQuantJob);

                    // new job
                    CondorJob sortBAMByReferenceAndNameJob = PipelineJobFactory.createJob(++count,
                            SortByReferenceAndNameCLI.class, getWorkflowPlan(), htsfSample);
                    sortBAMByReferenceAndNameJob.addArgument(SortByReferenceAndNameCLI.INPUT,
                            samtoolsSortOut.getAbsolutePath());
                    File sortBAMByReferenceAndNameOut = new File(outputDirectory, samtoolsSortOut.getName().replace(
                            ".bam", ".refAndNameSort.bam"));
                    sortBAMByReferenceAndNameJob.addArgument(SortByReferenceAndNameCLI.OUTPUT,
                            sortBAMByReferenceAndNameOut.getAbsolutePath());
                    sortBAMByReferenceAndNameJob
                            .addArgument(SortByReferenceAndNameCLI.TMPDIR, tmpDir.getAbsolutePath());
                    graph.addVertex(sortBAMByReferenceAndNameJob);
                    graph.addEdge(samtoolsIndexJob, sortBAMByReferenceAndNameJob);

                    // new job
                    CondorJob ubuSamTranslateJob = PipelineJobFactory.createJob(++count, UBUSamTranslateCLI.class,
                            getWorkflowPlan(), htsfSample);
                    ubuSamTranslateJob.addArgument(UBUSamTranslateCLI.BED, this.pipelineBeanService.getBed());
                    ubuSamTranslateJob.addArgument(UBUSamTranslateCLI.INPUT,
                            sortBAMByReferenceAndNameOut.getAbsolutePath());
                    File ubuSamTranslateOut = new File(outputDirectory, sortBAMByReferenceAndNameOut.getName().replace(
                            ".bam", ".transcriptomeAlignments.bam"));
                    ubuSamTranslateJob.addArgument(UBUSamTranslateCLI.OUTPUT, ubuSamTranslateOut.getAbsolutePath());
                    ubuSamTranslateJob.addArgument(UBUSamTranslateCLI.ORDER, this.pipelineBeanService.getOrder());
                    ubuSamTranslateJob.addArgument(UBUSamTranslateCLI.XGTAGS);
                    ubuSamTranslateJob.addArgument(UBUSamTranslateCLI.REVERSE);
                    graph.addVertex(ubuSamTranslateJob);
                    graph.addEdge(sortBAMByReferenceAndNameJob, ubuSamTranslateJob);

                    // new job
                    CondorJob ubuSamFilterJob = PipelineJobFactory.createJob(++count, UBUSamFilterCLI.class,
                            getWorkflowPlan(), htsfSample);
                    ubuSamFilterJob.addArgument(UBUSamFilterCLI.INPUT, ubuSamTranslateOut.getAbsolutePath());
                    File ubuSamFilterOut = new File(outputDirectory, ubuSamTranslateOut.getName().replace(".bam",
                            ".filtered.bam"));
                    ubuSamFilterJob.addArgument(UBUSamFilterCLI.OUTPUT, ubuSamFilterOut.getAbsolutePath());
                    ubuSamFilterJob.addArgument(UBUSamFilterCLI.STRIPINDELS);
                    ubuSamFilterJob.addArgument(UBUSamFilterCLI.MAXINSERT, "10000");
                    ubuSamFilterJob.addArgument(UBUSamFilterCLI.MAPQ, "1");
                    graph.addVertex(ubuSamFilterJob);
                    graph.addEdge(ubuSamTranslateJob, ubuSamFilterJob);

                    // new job
                    CondorJob rsemJob = PipelineJobFactory.createJob(++count, RSEMCalculateExpressionCLI.class,
                            getWorkflowPlan(), htsfSample);
                    rsemJob.addArgument(RSEMCalculateExpressionCLI.PAIREDEND);
                    rsemJob.addArgument(RSEMCalculateExpressionCLI.BAM);
                    rsemJob.addArgument(RSEMCalculateExpressionCLI.ESTIMATERSPD);
                    rsemJob.addArgument(RSEMCalculateExpressionCLI.THREADS, "4");
                    rsemJob.setNumberOfProcessors(4);
                    rsemJob.addArgument(RSEMCalculateExpressionCLI.BAMFILE, ubuSamFilterOut.getAbsolutePath());
                    rsemJob.addArgument(RSEMCalculateExpressionCLI.REFERENCESEQUENCE,
                            this.pipelineBeanService.getReferenceSequencePrefix());
                    File outputPrefix = new File(outputDirectory, "rsem");
                    rsemJob.addArgument(RSEMCalculateExpressionCLI.OUTPUT, outputPrefix.getAbsolutePath());
                    graph.addVertex(rsemJob);
                    graph.addEdge(ubuSamFilterJob, rsemJob);

                    // new job
                    CondorJob rsemISOFormsResultStripTabJob = PipelineJobFactory.createJob(++count,
                            StripTrailingTabsCLI.class, getWorkflowPlan(), htsfSample);
                    File isoformResults = new File(outputDirectory, outputPrefix.getName() + ".isoforms.results");
                    rsemISOFormsResultStripTabJob.addArgument(StripTrailingTabsCLI.INPUT,
                            isoformResults.getAbsolutePath());
                    File isoformResultsStripped = new File(outputDirectory, outputPrefix.getName()
                            + ".isoforms.stripped.results");
                    rsemISOFormsResultStripTabJob.addArgument(StripTrailingTabsCLI.OUTPUT,
                            isoformResultsStripped.getAbsolutePath());
                    rsemISOFormsResultStripTabJob.setPostScript(String.format("/bin/mv %s %s",
                            isoformResultsStripped.getAbsolutePath(), isoformResults.getAbsolutePath()));
                    graph.addVertex(rsemISOFormsResultStripTabJob);
                    graph.addEdge(rsemJob, rsemISOFormsResultStripTabJob);

                    // new job
                    CondorJob pruneISOFormsFromGeneQuantFileJob = PipelineJobFactory.createJob(++count,
                            PruneISOFormsFromGeneQuantFileCLI.class, getWorkflowPlan(), htsfSample);
                    File geneResultsFile = new File(outputDirectory, outputPrefix.getName() + ".genes.results");
                    pruneISOFormsFromGeneQuantFileJob.addArgument(PruneISOFormsFromGeneQuantFileCLI.GENERESULTS,
                            geneResultsFile.getAbsolutePath());
                    File origGeneResultsFile = new File(outputDirectory, "orig." + outputPrefix.getName()
                            + ".genes.results");
                    pruneISOFormsFromGeneQuantFileJob.addArgument(PruneISOFormsFromGeneQuantFileCLI.ORIGGENERESULTS,
                            origGeneResultsFile.getAbsolutePath());
                    graph.addVertex(pruneISOFormsFromGeneQuantFileJob);
                    graph.addEdge(rsemISOFormsResultStripTabJob, pruneISOFormsFromGeneQuantFileJob);

                    // new job
                    CondorJob normalizeGeneQuantJob = PipelineJobFactory.createJob(++count, NormalizeQuartileCLI.class,
                            getWorkflowPlan(), htsfSample);
                    normalizeGeneQuantJob.addArgument(NormalizeQuartileCLI.COLUMN, "2");
                    normalizeGeneQuantJob.addArgument(NormalizeQuartileCLI.QUANTILE, "75");
                    normalizeGeneQuantJob.addArgument(NormalizeQuartileCLI.TARGET, "1000");
                    normalizeGeneQuantJob.addArgument(NormalizeQuartileCLI.INPUT, geneResultsFile.getAbsolutePath());
                    File normalizeGeneQuantOut = new File(outputDirectory, outputPrefix.getName()
                            + ".genes.normalized_results");
                    normalizeGeneQuantJob.addArgument(NormalizeQuartileCLI.OUTPUT,
                            normalizeGeneQuantOut.getAbsolutePath());
                    graph.addVertex(normalizeGeneQuantJob);
                    graph.addEdge(pruneISOFormsFromGeneQuantFileJob, normalizeGeneQuantJob);

                    // new job
                    CondorJob normalizeISOFormQuantJob = PipelineJobFactory.createJob(++count,
                            NormalizeQuartileCLI.class, getWorkflowPlan(), htsfSample);
                    normalizeISOFormQuantJob.addArgument(NormalizeQuartileCLI.COLUMN, "2");
                    normalizeISOFormQuantJob.addArgument(NormalizeQuartileCLI.QUANTILE, "75");
                    normalizeISOFormQuantJob.addArgument(NormalizeQuartileCLI.TARGET, "300");
                    normalizeISOFormQuantJob.addArgument(NormalizeQuartileCLI.INPUT, isoformResults.getAbsolutePath());
                    File normalizeISOFormQuantOut = new File(outputDirectory, outputPrefix.getName()
                            + ".isoforms.normalized_results");
                    normalizeISOFormQuantJob.addArgument(NormalizeQuartileCLI.OUTPUT,
                            normalizeISOFormQuantOut.getAbsolutePath());
                    graph.addVertex(normalizeISOFormQuantJob);
                    graph.addEdge(pruneISOFormsFromGeneQuantFileJob, normalizeISOFormQuantJob);
                } catch (Exception e) {
                    throw new PipelineException(e);
                }

            }

        }

        return graph;
    }

    public RNASeqPipelineBeanService getPipelineBeanService() {
        return pipelineBeanService;
    }

    public void setPipelineBeanService(RNASeqPipelineBeanService pipelineBeanService) {
        this.pipelineBeanService = pipelineBeanService;
    }

}

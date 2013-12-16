package edu.unc.mapseq.workflow.rnaseq;

import java.io.File;
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

import edu.unc.mapseq.dao.model.HTSFSample;
import edu.unc.mapseq.dao.model.SequencerRun;
import edu.unc.mapseq.module.bedtools.CoverageBedCLI;
import edu.unc.mapseq.module.bowtie.BowtieBuildCLI;
import edu.unc.mapseq.module.core.CatCLI;
import edu.unc.mapseq.module.core.DetermineReadLengthCLI;
import edu.unc.mapseq.module.core.GUnZipCLI;
import edu.unc.mapseq.module.core.RegexCatCLI;
import edu.unc.mapseq.module.core.RemoveCLI;
import edu.unc.mapseq.module.core.SedCLI;
import edu.unc.mapseq.module.core.SortCLI;
import edu.unc.mapseq.module.filter.PruneISOFormsFromGeneQuantFileCLI;
import edu.unc.mapseq.module.filter.StripTrailingTabsCLI;
import edu.unc.mapseq.module.mapsplice.AlignmentHandlerMultiCLI;
import edu.unc.mapseq.module.mapsplice.ClusterCLI;
import edu.unc.mapseq.module.mapsplice.Filter1HitsCLI;
import edu.unc.mapseq.module.mapsplice.FilterJunctionByROCarguNonCanonicalCLI;
import edu.unc.mapseq.module.mapsplice.FilterOriginalFusionCLI;
import edu.unc.mapseq.module.mapsplice.FusionSAM2JunctionFilterAnchorNewFormatCLI;
import edu.unc.mapseq.module.mapsplice.JunctionDBFusionCLI;
import edu.unc.mapseq.module.mapsplice.JunctionSequenceConstructionCLI;
import edu.unc.mapseq.module.mapsplice.MapSpliceMultiThreadCLI;
import edu.unc.mapseq.module.mapsplice.NewSAM2JuntionCLI;
import edu.unc.mapseq.module.mapsplice.ParseClusterCLI;
import edu.unc.mapseq.module.mapsplice.QualityScaleType;
import edu.unc.mapseq.module.mapsplice.RSEMCalculateExpressionCLI;
import edu.unc.mapseq.module.mapsplice.ReadChromoSizeCLI;
import edu.unc.mapseq.module.mapsplice.ReadsToUnmappedSAMCLI;
import edu.unc.mapseq.module.mapsplice.SetUnmappedBitFlagCLI;
import edu.unc.mapseq.module.picard.PicardAddOrReplaceReadGroupsCLI;
import edu.unc.mapseq.module.picard.PicardSortOrderType;
import edu.unc.mapseq.module.qc.NormBedExonQuantCLI;
import edu.unc.mapseq.module.qc.NormalizeQuartileCLI;
import edu.unc.mapseq.module.samtools.SAMToolsFlagstatCLI;
import edu.unc.mapseq.module.samtools.SAMToolsIndexCLI;
import edu.unc.mapseq.module.samtools.SAMToolsSortCLI;
import edu.unc.mapseq.module.samtools.SAMToolsViewCLI;
import edu.unc.mapseq.module.samtools.SortByReferenceAndNameCLI;
import edu.unc.mapseq.module.ubu.UBUFastqFormatterCLI;
import edu.unc.mapseq.module.ubu.UBUSamFilterCLI;
import edu.unc.mapseq.module.ubu.UBUSamJunctionCLI;
import edu.unc.mapseq.module.ubu.UBUSamTranslateCLI;
import edu.unc.mapseq.workflow.AbstractWorkflow;
import edu.unc.mapseq.workflow.WorkflowException;
import edu.unc.mapseq.workflow.WorkflowJobFactory;
import edu.unc.mapseq.workflow.WorkflowUtil;

public class RNASeqWorkflow extends AbstractWorkflow {

    private final Logger logger = LoggerFactory.getLogger(RNASeqWorkflow.class);

    public RNASeqWorkflow() {
        super();
    }

    @Override
    public String getName() {
        return RNASeqWorkflow.class.getSimpleName().replace("Workflow", "");
    }

    @Override
    public String getVersion() {
        ResourceBundle bundle = ResourceBundle.getBundle("edu/unc/mapseq/workflow/rnaseq/workflow");
        String version = bundle.getString("version");
        return StringUtils.isNotEmpty(version) ? version : "0.0.1-SNAPSHOT";
    }

    @Override
    public void preRun() throws WorkflowException {

        Set<HTSFSample> htsfSampleSet = getAggregateHTSFSampleSet();
        logger.info("htsfSampleSet.size(): {}", htsfSampleSet.size());

        for (HTSFSample htsfSample : htsfSampleSet) {

            if ("Undetermined".equals(htsfSample.getBarcode())) {
                continue;
            }

            SequencerRun sequencerRun = htsfSample.getSequencerRun();
            File outputDirectory = createOutputDirectory(sequencerRun.getName(), htsfSample, getName(), getVersion());

            File tmpDir = new File(outputDirectory, "tmp");

            File originalDir = new File(tmpDir, "original");
            originalDir.mkdirs();

            File bestDir = new File(tmpDir, "best");
            bestDir.mkdirs();

            File remapDir = new File(tmpDir, "remap");
            remapDir.mkdirs();

            File fusionDir = new File(tmpDir, "fusion");

            File fusionDataDir = new File(fusionDir, "data");

            File fusionDataPERDir = new File(fusionDataDir, "PER");
            fusionDataPERDir.mkdirs();

            File fusionDataSingleDir = new File(fusionDataDir, "single");
            fusionDataSingleDir.mkdirs();

            File fusionResultDir = new File(fusionDir, "result");

            File fusionResultPERProbDir = new File(fusionResultDir, "PER_prob");
            fusionResultPERProbDir.mkdirs();

            File fusionResultFusionReadDir = new File(fusionResultDir, "fusionRead");
            fusionResultFusionReadDir.mkdirs();

            File fusionResultJunctionSupportDir = new File(fusionResultDir, "junction_support");
            fusionResultJunctionSupportDir.mkdirs();

            File clusterDir = new File(tmpDir, "cluster");
            clusterDir.mkdirs();

            File resultDir = new File(clusterDir, "result");

            File pairedEndMatchDir = new File(resultDir, "PEmatch");
            pairedEndMatchDir.mkdirs();

            File pairedEndMatchRegionBEDDir = new File(resultDir, "PEmatch_region_BED");
            pairedEndMatchRegionBEDDir.mkdirs();

            File pairedEndMatchRegionDir = new File(resultDir, "PEmatchregion");
            pairedEndMatchRegionDir.mkdirs();

            File pairedEndMatchRegionDetailDir = new File(resultDir, "PEmatchregion_detail");
            pairedEndMatchRegionDetailDir.mkdirs();

            File distanceDir = new File(resultDir, "distance");
            distanceDir.mkdirs();

        }

    }

    @Override
    public Graph<CondorJob, CondorJobEdge> createGraph() throws WorkflowException {
        logger.info("ENTERING createGraph()");

        DirectedGraph<CondorJob, CondorJobEdge> graph = new DefaultDirectedGraph<CondorJob, CondorJobEdge>(
                CondorJobEdge.class);

        int count = 0;

        Set<HTSFSample> htsfSampleSet = getAggregateHTSFSampleSet();
        logger.info("htsfSampleSet.size(): {}", htsfSampleSet.size());

        String siteName = getWorkflowBeanService().getAttributes().get("siteName");

        String bowtieIndexDirectory = getWorkflowBeanService().getAttributes().get("bowtieIndexDirectory");
        logger.info("bowtieIndexDirectory: {}", bowtieIndexDirectory);

        String chromosomeDirectory = getWorkflowBeanService().getAttributes().get("chromosomeDirectory");
        logger.info("chromosomeDirectory: {}", chromosomeDirectory);

        String junctions = getWorkflowBeanService().getAttributes().get("junctions");
        logger.info("junctions: {}", junctions);

        String compositeExons = getWorkflowBeanService().getAttributes().get("compositeExons");
        logger.info("compositeExons: {}", compositeExons);

        String bed = getWorkflowBeanService().getAttributes().get("bed");
        logger.info("bed: {}", bed);

        String referenceSequencePrefix = getWorkflowBeanService().getAttributes().get("referenceSequencePrefix");
        logger.info("referenceSequencePrefix: {}", referenceSequencePrefix);

        for (HTSFSample htsfSample : htsfSampleSet) {

            if ("Undetermined".equals(htsfSample.getBarcode())) {
                continue;
            }

            logger.debug("htsfSample = {}", htsfSample.toString());

            SequencerRun sequencerRun = htsfSample.getSequencerRun();
            File outputDirectory = createOutputDirectory(sequencerRun.getName(), htsfSample, getName(), getVersion());

            File tmpDir = new File(outputDirectory, "tmp");
            File originalDir = new File(tmpDir, "original");
            File bestDir = new File(tmpDir, "best");
            File remapDir = new File(tmpDir, "remap");
            File fusionDir = new File(tmpDir, "fusion");
            File clusterDir = new File(tmpDir, "cluster");
            File resultDir = new File(clusterDir, "result");

            List<File> readPairList = WorkflowUtil.getReadPairList(htsfSample.getFileDatas(), sequencerRun.getName(),
                    htsfSample.getLaneIndex());
            logger.info("fileList = {}", readPairList.size());

            // assumption: a dash is used as a delimiter between a participantId
            // and the external code
            int idx = htsfSample.getName().lastIndexOf("-");
            String sampleName = idx != -1 ? htsfSample.getName().substring(0, idx) : htsfSample.getName();

            if (readPairList.size() != 2) {
                logger.error("readPairList.size() != 2");
                throw new WorkflowException("Need two fastq files to process...check FileData/Entity mapping");
            }

            File r1FastqFile = readPairList.get(0);
            String r1FastqRootName = WorkflowUtil.getRootFastqName(r1FastqFile.getName());

            File r2FastqFile = readPairList.get(1);
            String r2FastqRootName = WorkflowUtil.getRootFastqName(r2FastqFile.getName());

            String fastqLaneRootName = StringUtils.removeEnd(r2FastqRootName, "_R2");

            // new job
            CondorJob gunzipFastqR1Job = WorkflowJobFactory.createJob(++count, GUnZipCLI.class, getWorkflowPlan(),
                    htsfSample);
            gunzipFastqR1Job.setSiteName(siteName);
            gunzipFastqR1Job.addArgument(GUnZipCLI.GZFILE, r1FastqFile.getAbsolutePath());
            File gunzippedFastqR1 = new File(outputDirectory, r1FastqRootName + ".fastq");
            gunzipFastqR1Job.addArgument(GUnZipCLI.EXTRACTFILE, gunzippedFastqR1.getAbsolutePath());
            graph.addVertex(gunzipFastqR1Job);

            // new job
            CondorJob fastqFormatterR1Job = WorkflowJobFactory.createJob(++count, UBUFastqFormatterCLI.class,
                    getWorkflowPlan(), htsfSample);
            fastqFormatterR1Job.setSiteName(siteName);
            fastqFormatterR1Job.addArgument(UBUFastqFormatterCLI.INPUT, gunzippedFastqR1.getAbsolutePath());
            fastqFormatterR1Job.addArgument(UBUFastqFormatterCLI.STRIP);
            fastqFormatterR1Job.addArgument(UBUFastqFormatterCLI.SUFFIX, "/1");
            File fastqFormatterR1Out = new File(outputDirectory, gunzippedFastqR1.getName().replace(".fastq",
                    ".filtered.fastq"));
            fastqFormatterR1Job.addArgument(UBUFastqFormatterCLI.OUTPUT, fastqFormatterR1Out.getAbsolutePath());
            graph.addVertex(fastqFormatterR1Job);
            graph.addEdge(gunzipFastqR1Job, fastqFormatterR1Job);

            // new job
            CondorJob gunzipFastqR2Job = WorkflowJobFactory.createJob(++count, GUnZipCLI.class, getWorkflowPlan(),
                    htsfSample);
            gunzipFastqR2Job.setSiteName(siteName);
            gunzipFastqR2Job.addArgument(GUnZipCLI.GZFILE, r2FastqFile.getAbsolutePath());
            File gunzippedFastqR2 = new File(outputDirectory, r2FastqRootName + ".fastq");
            gunzipFastqR2Job.addArgument(GUnZipCLI.EXTRACTFILE, gunzippedFastqR2.getAbsolutePath());
            graph.addVertex(gunzipFastqR2Job);

            // new job
            CondorJob fastqFormatterR2Job = WorkflowJobFactory.createJob(++count, UBUFastqFormatterCLI.class,
                    getWorkflowPlan(), htsfSample);
            fastqFormatterR2Job.setSiteName(siteName);
            fastqFormatterR2Job.addArgument(UBUFastqFormatterCLI.INPUT, gunzippedFastqR2.getAbsolutePath());
            fastqFormatterR2Job.addArgument(UBUFastqFormatterCLI.STRIP);
            fastqFormatterR2Job.addArgument(UBUFastqFormatterCLI.SUFFIX, "/2");
            File fastqFormatterR2Out = new File(outputDirectory, gunzippedFastqR2.getName().replace(".fastq",
                    ".filtered.fastq"));
            fastqFormatterR2Job.addArgument(UBUFastqFormatterCLI.OUTPUT, fastqFormatterR2Out.getAbsolutePath());
            graph.addVertex(fastqFormatterR2Job);
            graph.addEdge(gunzipFastqR2Job, fastqFormatterR2Job);

            // new job
            CondorJob determineReadLengthJob = WorkflowJobFactory.createJob(++count, DetermineReadLengthCLI.class,
                    getWorkflowPlan(), htsfSample);
            determineReadLengthJob.setSiteName(siteName);
            determineReadLengthJob.addArgument(DetermineReadLengthCLI.INPUT, fastqFormatterR1Out.getAbsolutePath());
            determineReadLengthJob.addArgument(DetermineReadLengthCLI.INPUT, fastqFormatterR2Out.getAbsolutePath());
            File determineReadLengthOutput = new File(tmpDir, "readLengthProps.xml");
            determineReadLengthJob.addArgument(DetermineReadLengthCLI.OUTPUT,
                    determineReadLengthOutput.getAbsolutePath());
            graph.addVertex(determineReadLengthJob);
            graph.addEdge(fastqFormatterR1Job, determineReadLengthJob);
            graph.addEdge(fastqFormatterR2Job, determineReadLengthJob);

            // new job
            CondorJob readChromoSizesJob = WorkflowJobFactory.createJob(++count, ReadChromoSizeCLI.class,
                    getWorkflowPlan(), htsfSample);
            readChromoSizesJob.setSiteName(siteName);
            File chromosomeHeadOutput = new File(tmpDir, "chrom_head");
            readChromoSizesJob.addArgument(ReadChromoSizeCLI.CHROMOSOMEHEADOUTPUT,
                    chromosomeHeadOutput.getAbsolutePath());
            File chromosomeIndexOutput = new File(tmpDir, "chromo.fai");
            readChromoSizesJob.addArgument(ReadChromoSizeCLI.CHROMOSOMEINDEX, chromosomeIndexOutput.getAbsolutePath());
            File chromosomeNamesOutput = new File(fusionDir, "chrName.txt");
            readChromoSizesJob.addArgument(ReadChromoSizeCLI.CHROMOSOMENAMESOUTPUT,
                    chromosomeNamesOutput.getAbsolutePath());
            File chromosomeSizesOutput = new File(tmpDir, "chrom_sizes");
            readChromoSizesJob.addArgument(ReadChromoSizeCLI.CHROMOSOMESIZESOUTPUT,
                    chromosomeSizesOutput.getAbsolutePath());
            File chromDir = new File(chromosomeDirectory);
            for (File f : chromDir.listFiles()) {
                if (f.getName().endsWith(".fa")) {
                    readChromoSizesJob.addArgument(ReadChromoSizeCLI.SAMFILE, f.getAbsolutePath());
                }
            }
            graph.addVertex(readChromoSizesJob);

            // new job
            CondorJob originalSAMMapspliceMultiThreadJob = WorkflowJobFactory.createJob(++count,
                    MapSpliceMultiThreadCLI.class, getWorkflowPlan(), htsfSample);
            originalSAMMapspliceMultiThreadJob.setSiteName(siteName);
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.FASTQ1,
                    fastqFormatterR1Out.getAbsolutePath());
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.FASTQ2,
                    fastqFormatterR2Out.getAbsolutePath());
            originalSAMMapspliceMultiThreadJob.setNumberOfProcessors(8);
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.THREADS, 8);
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.QUALITYSCALE,
                    QualityScaleType.phred33);
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.MINIMUMINTRON, 50);
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.SUPPRESSALIGNMENTSOVER, 40);
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMINTRONSINGLE, 200000);
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMINTRONDOUBLE, 200000);
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.MINIMUMLENGTH, 25);
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.SPLICEMISMATCHES, 1);
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.SPLICEONLY);
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.REPORTUPTONALIGNMENTSPERREAD, 40);
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.FASTQREADFORMAT);
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMINSERTIONS, 6);
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMDELETIONS, 6);
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMAPPENDMISMATCHES, 3);
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.SEGMENTLENGTH, 25);
            File mapspliceMultiThreadJobOutput = new File(originalDir, "bowtie.output");
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.OUTPUT,
                    mapspliceMultiThreadJobOutput.getAbsolutePath());
            File originalSAMMapspliceMultiThreadOutput = new File(originalDir, "original.sam");
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.MAPSPLICEOUT,
                    originalSAMMapspliceMultiThreadOutput.getAbsolutePath());
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.REFERENCESEQUENCEDIRECTORY,
                    chromosomeDirectory);
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.CHROMOSOMESIZE,
                    chromosomeSizesOutput.getAbsolutePath());
            originalSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.INDEX, bowtieIndexDirectory);
            graph.addVertex(originalSAMMapspliceMultiThreadJob);
            graph.addEdge(fastqFormatterR1Job, originalSAMMapspliceMultiThreadJob);
            graph.addEdge(fastqFormatterR2Job, originalSAMMapspliceMultiThreadJob);

            // new job
            CondorJob originalSAMRegexCatJob = WorkflowJobFactory.createJob(++count, RegexCatCLI.class,
                    getWorkflowPlan(), htsfSample, false);
            originalSAMRegexCatJob.setSiteName(siteName);
            originalSAMRegexCatJob.addArgument(RegexCatCLI.DIRECTORY, originalDir.getAbsolutePath());
            originalSAMRegexCatJob.addArgument(RegexCatCLI.OUTPUT,
                    originalSAMMapspliceMultiThreadOutput.getAbsolutePath());
            originalSAMRegexCatJob.addArgument(RegexCatCLI.REGEX, "^.*\\.sam\\.[0-9]");
            graph.addVertex(originalSAMRegexCatJob);
            graph.addEdge(originalSAMMapspliceMultiThreadJob, originalSAMRegexCatJob);

            // new job
            CondorJob sam2JunctionArrayJob = WorkflowJobFactory.createJob(++count, NewSAM2JuntionCLI.class,
                    getWorkflowPlan(), htsfSample, false);
            sam2JunctionArrayJob.setSiteName(siteName);
            sam2JunctionArrayJob.addArgument(NewSAM2JuntionCLI.CHROMOSOMEFILESDIRECTORY, chromosomeDirectory);
            sam2JunctionArrayJob.addArgument(NewSAM2JuntionCLI.INPUT,
                    originalSAMMapspliceMultiThreadOutput.getAbsolutePath());
            File originalAllJunctionsFile = new File(originalDir, "ori.all_junctions.txt");
            sam2JunctionArrayJob
                    .addArgument(NewSAM2JuntionCLI.JUNCTIONFILE, originalAllJunctionsFile.getAbsolutePath());
            sam2JunctionArrayJob.addArgument(NewSAM2JuntionCLI.MINIMUMANCHOR, 1);
            sam2JunctionArrayJob.addArgument(NewSAM2JuntionCLI.MINIMUMINTRON, 1);
            sam2JunctionArrayJob.addArgument(NewSAM2JuntionCLI.MAXIMUMINTRON, 350000000);
            sam2JunctionArrayJob.addArgument(NewSAM2JuntionCLI.READLENGTH, determineReadLengthOutput.getAbsolutePath());
            graph.addVertex(sam2JunctionArrayJob);
            graph.addEdge(determineReadLengthJob, sam2JunctionArrayJob);
            graph.addEdge(originalSAMRegexCatJob, sam2JunctionArrayJob);

            // new job
            CondorJob filter1HitsJob = WorkflowJobFactory.createJob(++count, Filter1HitsCLI.class, getWorkflowPlan(),
                    htsfSample, false);
            filter1HitsJob.setSiteName(siteName);
            filter1HitsJob.addArgument(Filter1HitsCLI.JUNCTIONFILE, originalAllJunctionsFile.getAbsolutePath());
            File canonicalFilter1HitsOutput = new File(bestDir,
                    "ori.all_junctions.filtered_by_min_mis_lpq.remained.txt");
            filter1HitsJob.addArgument(Filter1HitsCLI.CANONICALJUNCTIONOUTPUT,
                    canonicalFilter1HitsOutput.getAbsolutePath());
            File nonCanonicalFilter1HitsOutput = new File(bestDir,
                    "ori.all_junctions.filtered_by_min_mis_lpq.filtered.txt");
            filter1HitsJob.addArgument(Filter1HitsCLI.NONCANONICALJUNCTIONOUTPUT,
                    nonCanonicalFilter1HitsOutput.getAbsolutePath());
            filter1HitsJob.addArgument(Filter1HitsCLI.MINIMUMLPQ, 0.3);
            filter1HitsJob.addArgument(Filter1HitsCLI.MINIMUMMISMATCH, 2);
            graph.addVertex(filter1HitsJob);
            graph.addEdge(sam2JunctionArrayJob, filter1HitsJob);

            // new job
            CondorJob filterJunctionByROCarguNonCanonicalJob = WorkflowJobFactory.createJob(++count,
                    FilterJunctionByROCarguNonCanonicalCLI.class, getWorkflowPlan(), htsfSample, false);
            filterJunctionByROCarguNonCanonicalJob.setSiteName(siteName);
            filterJunctionByROCarguNonCanonicalJob.addArgument(FilterJunctionByROCarguNonCanonicalCLI.ENTROPYWEIGHT,
                    0.097718);
            filterJunctionByROCarguNonCanonicalJob.addArgument(FilterJunctionByROCarguNonCanonicalCLI.LPQWEIGHT,
                    0.66478);
            filterJunctionByROCarguNonCanonicalJob.addArgument(FilterJunctionByROCarguNonCanonicalCLI.MINSCORE, 0.719);
            filterJunctionByROCarguNonCanonicalJob.addArgument(
                    FilterJunctionByROCarguNonCanonicalCLI.AVEMISMATCHWEIGHT, -0.21077);
            filterJunctionByROCarguNonCanonicalJob.addArgument(FilterJunctionByROCarguNonCanonicalCLI.JUNCTIONFILE,
                    canonicalFilter1HitsOutput.getAbsolutePath());
            File bestJunctionOutput = new File(bestDir, "best_junction.txt");
            filterJunctionByROCarguNonCanonicalJob.addArgument(
                    FilterJunctionByROCarguNonCanonicalCLI.CANONICALJUNCTIONOUTPUT,
                    bestJunctionOutput.getAbsolutePath());
            File nonCanonicalBestJunctionOutput = new File(bestDir,
                    "best_junction_semi_non_canon_filtered_by_ROCargu.txt");
            filterJunctionByROCarguNonCanonicalJob.addArgument(
                    FilterJunctionByROCarguNonCanonicalCLI.NONCANONICALJUNCTIONOUTPUT,
                    nonCanonicalBestJunctionOutput.getAbsolutePath());
            filterJunctionByROCarguNonCanonicalJob.addArgument(FilterJunctionByROCarguNonCanonicalCLI.INTRONWEIGHT, 0);
            filterJunctionByROCarguNonCanonicalJob.addArgument(FilterJunctionByROCarguNonCanonicalCLI.SUMLENGTHWEIGHT,
                    0);
            filterJunctionByROCarguNonCanonicalJob.addArgument(FilterJunctionByROCarguNonCanonicalCLI.MINFLANKCASE, 5);
            graph.addVertex(filterJunctionByROCarguNonCanonicalJob);
            graph.addEdge(filter1HitsJob, filterJunctionByROCarguNonCanonicalJob);

            // new job
            CondorJob junctionSequenceConstructionJob = WorkflowJobFactory.createJob(++count,
                    JunctionSequenceConstructionCLI.class, getWorkflowPlan(), htsfSample, false);
            junctionSequenceConstructionJob.setSiteName(siteName);
            junctionSequenceConstructionJob.addArgument(JunctionSequenceConstructionCLI.JUNCTION,
                    bestJunctionOutput.getAbsolutePath());
            junctionSequenceConstructionJob.addArgument(JunctionSequenceConstructionCLI.REFERENCESEQUENCEDIRECTORY,
                    chromosomeDirectory);
            File junctionSequenceConstructionOutput = new File(remapDir, "synthetic_alljunc_sequence.txt");
            junctionSequenceConstructionJob.addArgument(JunctionSequenceConstructionCLI.OUTPUT,
                    junctionSequenceConstructionOutput.getAbsolutePath());
            junctionSequenceConstructionJob.addArgument(JunctionSequenceConstructionCLI.MINIMUMANCHOR, 2);
            junctionSequenceConstructionJob.addArgument(JunctionSequenceConstructionCLI.MAXIMUMANCHOR, 38);
            junctionSequenceConstructionJob.addArgument(JunctionSequenceConstructionCLI.MAXIMUMSEQUENCETHRESHOLD, 400);
            graph.addVertex(junctionSequenceConstructionJob);
            graph.addEdge(filterJunctionByROCarguNonCanonicalJob, junctionSequenceConstructionJob);

            CondorJob bowtieBuildJob = WorkflowJobFactory.createJob(++count, BowtieBuildCLI.class, getWorkflowPlan(),
                    htsfSample, false);
            bowtieBuildJob.setSiteName(siteName);
            bowtieBuildJob.addArgument(BowtieBuildCLI.INPUT, junctionSequenceConstructionOutput.getAbsolutePath());
            File synIndexPrefix = new File(remapDir, "syn_idx_prefix");
            bowtieBuildJob.addArgument(BowtieBuildCLI.PREFIX, synIndexPrefix.getAbsolutePath());
            graph.addVertex(bowtieBuildJob);
            graph.addEdge(junctionSequenceConstructionJob, bowtieBuildJob);

            // new job
            CondorJob remappedSAMMapspliceMultiThreadJob = WorkflowJobFactory.createJob(++count,
                    MapSpliceMultiThreadCLI.class, getWorkflowPlan(), htsfSample);
            remappedSAMMapspliceMultiThreadJob.setSiteName(siteName);
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.FASTQ1,
                    fastqFormatterR1Out.getAbsolutePath());
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.FASTQ2,
                    fastqFormatterR2Out.getAbsolutePath());
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.THREADS, 8);
            remappedSAMMapspliceMultiThreadJob.setNumberOfProcessors(8);
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.QUALITYSCALE,
                    QualityScaleType.phred33);
            File unmappedOutput = new File(remapDir, "remap_unmapped");
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.OUTPUTUNMAPPED,
                    unmappedOutput.getAbsolutePath());
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.MINIMUMINTRON, 50);
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMINTRONSINGLE, 200000);
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMINTRONDOUBLE, 200000);
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.MINIMUMLENGTH, 25);
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.SEGMENTLENGTH, 25);
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.SUPPRESSALIGNMENTSOVER, 40);
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.REPORTUPTONALIGNMENTSPERREAD, 40);
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.OPTIMIZEREPEATS);
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.JUNCTIONINDEX,
                    synIndexPrefix.getAbsolutePath());
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.SPLICEMISMATCHES, 1);
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.FASTQREADFORMAT);
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMINSERTIONS, 6);
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMDELETIONS, 6);
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMAPPENDMISMATCHES, 3);
            mapspliceMultiThreadJobOutput = new File(remapDir, "bowtie.output");
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.OUTPUT,
                    mapspliceMultiThreadJobOutput.getAbsolutePath());
            File remappedSAMMapspliceMultiThreadOutput = new File(remapDir, "remapped.sam");
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.MAPSPLICEOUT,
                    remappedSAMMapspliceMultiThreadOutput.getAbsolutePath());
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.REFERENCESEQUENCEDIRECTORY,
                    chromosomeDirectory);
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.CHROMOSOMESIZE,
                    chromosomeSizesOutput.getAbsolutePath());
            remappedSAMMapspliceMultiThreadJob.addArgument(MapSpliceMultiThreadCLI.INDEX, bowtieIndexDirectory);
            graph.addVertex(remappedSAMMapspliceMultiThreadJob);
            graph.addEdge(bowtieBuildJob, remappedSAMMapspliceMultiThreadJob);

            // new job
            CondorJob remappedRegexCatJob = WorkflowJobFactory.createJob(++count, RegexCatCLI.class, getWorkflowPlan(),
                    htsfSample, false);
            remappedRegexCatJob.setSiteName(siteName);
            remappedRegexCatJob.addArgument(RegexCatCLI.DIRECTORY, remapDir.getAbsolutePath());
            remappedRegexCatJob
                    .addArgument(RegexCatCLI.OUTPUT, remappedSAMMapspliceMultiThreadOutput.getAbsolutePath());
            remappedRegexCatJob.addArgument(RegexCatCLI.REGEX, "^remapped\\.sam\\.[0-9]");
            graph.addVertex(remappedRegexCatJob);
            graph.addEdge(remappedSAMMapspliceMultiThreadJob, remappedRegexCatJob);

            // new job
            CondorJob remapUnmapped1RegexCatJob = WorkflowJobFactory.createJob(++count, RegexCatCLI.class,
                    getWorkflowPlan(), htsfSample, false);
            remapUnmapped1RegexCatJob.setSiteName(siteName);
            remapUnmapped1RegexCatJob.addArgument(RegexCatCLI.DIRECTORY, remapDir.getAbsolutePath());
            File remapUnmapped1RegexCatOutput = new File(remapDir, "remap_unmapped.1");
            remapUnmapped1RegexCatJob.addArgument(RegexCatCLI.OUTPUT, remapUnmapped1RegexCatOutput.getAbsolutePath());
            remapUnmapped1RegexCatJob.addArgument(RegexCatCLI.REGEX, "^remap_unmapped\\.[0-9]\\.1");
            graph.addVertex(remapUnmapped1RegexCatJob);
            graph.addEdge(remappedSAMMapspliceMultiThreadJob, remapUnmapped1RegexCatJob);

            // new job
            CondorJob remapUnmapped2RegexCatJob = WorkflowJobFactory.createJob(++count, RegexCatCLI.class,
                    getWorkflowPlan(), htsfSample, false);
            remapUnmapped2RegexCatJob.setSiteName(siteName);
            remapUnmapped2RegexCatJob.addArgument(RegexCatCLI.DIRECTORY, remapDir.getAbsolutePath());
            File remapUnmapped2RegexCatOutput = new File(remapDir, "remap_unmapped.2");
            remapUnmapped2RegexCatJob.addArgument(RegexCatCLI.OUTPUT, remapUnmapped2RegexCatOutput.getAbsolutePath());
            remapUnmapped2RegexCatJob.addArgument(RegexCatCLI.REGEX, "^remap_unmapped\\.[0-9]\\.2");
            graph.addVertex(remapUnmapped2RegexCatJob);
            graph.addEdge(remappedSAMMapspliceMultiThreadJob, remapUnmapped2RegexCatJob);

            // new job
            CondorJob alignmentHandlerMultiJob = WorkflowJobFactory.createJob(++count, AlignmentHandlerMultiCLI.class,
                    getWorkflowPlan(), htsfSample, false);
            alignmentHandlerMultiJob.setSiteName(siteName);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.ADDSOFTCLIP, 1);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.AVERAGEFRAGMENTLENGTH, 225);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.BOUNDARY, 36);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.CHROMOSOMESIZEFILE,
                    chromosomeSizesOutput.getAbsolutePath());
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.ENCOMPASSINGFUSIONREGIONEXTENSION, 50000);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.FILTEREDALIGNMENTAPPEND,
                    "_filtered_normal_alignments");
            File filteredAlignmentBase = new File(remapDir, "_filtered_normal_alignments");
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.FILTEREDALIGNMENTBASE,
                    filteredAlignmentBase.getAbsolutePath());
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.FILTERFLAG, 12 + 32 + 256);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.FRAGMENTLENGTH, 400);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.FRAGMENTLENGTHSD, 100);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.INPUT,
                    remappedSAMMapspliceMultiThreadOutput.getAbsolutePath());
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.INTRONDISTANCESD, 500);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MATEDISTANCESD, 100);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MAXIMUMANCHORDIFF, 50);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MAXIMUMDELETION, 6);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MAXIMUMHITS, 40);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MAXIMUMMATEDISTANCE, 50000);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.READLENGTHPROPERTIES,
                    determineReadLengthOutput.getAbsolutePath());
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MINIMUMANCHOR, 0);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MINIMUMCOVERAGE, 0);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MINIMUMENCOMPASSCOUNT, 1);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MINIMUMENTROPY, -0.0001);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MINIMUMINSERTION, 6);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MINIMUMJUNCTIONANCHOR, 10);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MINIMUMMISMATCH, 5);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.THREADS, 8);
            alignmentHandlerMultiJob.setNumberOfProcessors(8);
            graph.addVertex(alignmentHandlerMultiJob);
            graph.addEdge(determineReadLengthJob, alignmentHandlerMultiJob);
            graph.addEdge(remappedRegexCatJob, alignmentHandlerMultiJob);

            // new job
            CondorJob filteredNormalAlignmentsFusionPairedCatJob = WorkflowJobFactory.createJob(++count, CatCLI.class,
                    getWorkflowPlan(), htsfSample, false);
            filteredNormalAlignmentsFusionPairedCatJob.setSiteName(siteName);
            File fusionPaired = new File(remapDir, remappedSAMMapspliceMultiThreadOutput.getName()
                    + "_filtered_normal_alignments.fusion_paired");
            filteredNormalAlignmentsFusionPairedCatJob.addArgument(CatCLI.FILES, fusionPaired.getAbsolutePath());
            File bothUnsplicedFusionPairedAlignments = new File(remapDir,
                    remappedSAMMapspliceMultiThreadOutput.getName()
                            + "_filtered_normal_alignments.bothunspliced.fusion_paired");
            filteredNormalAlignmentsFusionPairedCatJob.addArgument(CatCLI.FILES,
                    bothUnsplicedFusionPairedAlignments.getAbsolutePath());
            File fusionPairedAlignments = new File(remapDir, remappedSAMMapspliceMultiThreadOutput.getName()
                    + "_filtered_normal_alignments.fusion_paired.comb");
            filteredNormalAlignmentsFusionPairedCatJob.addArgument(CatCLI.OUTPUT,
                    fusionPairedAlignments.getAbsolutePath());
            graph.addVertex(filteredNormalAlignmentsFusionPairedCatJob);
            graph.addEdge(alignmentHandlerMultiJob, filteredNormalAlignmentsFusionPairedCatJob);

            // new job
            CondorJob filteredNormalAlignmentsSingleCatJob = WorkflowJobFactory.createJob(++count, CatCLI.class,
                    getWorkflowPlan(), htsfSample, false);
            filteredNormalAlignmentsSingleCatJob.setSiteName(siteName);
            File single = new File(remapDir, remappedSAMMapspliceMultiThreadOutput.getName()
                    + "_filtered_normal_alignments.single");
            filteredNormalAlignmentsSingleCatJob.addArgument(CatCLI.FILES, single.getAbsolutePath());
            File bothUnsplicedSingleAlignments = new File(remapDir, remappedSAMMapspliceMultiThreadOutput.getName()
                    + "_filtered_normal_alignments.bothunspliced.single");
            filteredNormalAlignmentsSingleCatJob.addArgument(CatCLI.FILES,
                    bothUnsplicedSingleAlignments.getAbsolutePath());
            File singleAlignments = new File(remapDir, remappedSAMMapspliceMultiThreadOutput.getName()
                    + "_filtered_normal_alignments.single.comb");
            filteredNormalAlignmentsSingleCatJob.addArgument(CatCLI.OUTPUT, singleAlignments.getAbsolutePath());
            graph.addVertex(filteredNormalAlignmentsSingleCatJob);
            graph.addEdge(alignmentHandlerMultiJob, filteredNormalAlignmentsSingleCatJob);

            // new job
            CondorJob filteredNormalAlignmentsPairedCatJob = WorkflowJobFactory.createJob(++count, CatCLI.class,
                    getWorkflowPlan(), htsfSample, false);
            filteredNormalAlignmentsPairedCatJob.setSiteName(siteName);
            File paired = new File(remapDir, remappedSAMMapspliceMultiThreadOutput.getName()
                    + "_filtered_normal_alignments.paired");
            filteredNormalAlignmentsPairedCatJob.addArgument(CatCLI.FILES, paired.getAbsolutePath());
            File bothUnsplicedPairedAlignments = new File(remapDir, remappedSAMMapspliceMultiThreadOutput.getName()
                    + "_filtered_normal_alignments.bothunspliced.paired");
            filteredNormalAlignmentsPairedCatJob.addArgument(CatCLI.FILES,
                    bothUnsplicedPairedAlignments.getAbsolutePath());
            File pairedAlignments = new File(remapDir, remappedSAMMapspliceMultiThreadOutput.getName()
                    + "_filtered_normal_alignments.paired.comb");
            filteredNormalAlignmentsPairedCatJob.addArgument(CatCLI.OUTPUT, pairedAlignments.getAbsolutePath());
            graph.addVertex(filteredNormalAlignmentsPairedCatJob);
            graph.addEdge(alignmentHandlerMultiJob, filteredNormalAlignmentsPairedCatJob);

            // new job
            CondorJob sedJob = WorkflowJobFactory
                    .createJob(++count, SedCLI.class, getWorkflowPlan(), htsfSample, false);
            sedJob.setSiteName(siteName);
            sedJob.addArgument(SedCLI.SOURCE, fusionPairedAlignments.getAbsolutePath());
            sedJob.addArgument(SedCLI.REGULAREXPRESSION, "s/^/1~/");
            File fusionSAM = new File(fusionDir, "sed.fusion.sam");
            sedJob.addArgument(SedCLI.OUTPUT, fusionSAM.getAbsolutePath());
            graph.addVertex(sedJob);
            graph.addEdge(filteredNormalAlignmentsFusionPairedCatJob, sedJob);

            // new job
            CondorJob parseClusterJob = WorkflowJobFactory.createJob(++count, ParseClusterCLI.class, getWorkflowPlan(),
                    htsfSample, false);
            parseClusterJob.setSiteName(siteName);
            parseClusterJob.addArgument(ParseClusterCLI.INPUT, fusionSAM.getAbsolutePath());
            parseClusterJob.addArgument(ParseClusterCLI.OUTPUTDIRECTORY, clusterDir.getAbsolutePath());
            graph.addVertex(parseClusterJob);
            graph.addEdge(sedJob, parseClusterJob);

            // new job
            CondorJob clusterJob = WorkflowJobFactory.createJob(++count, ClusterCLI.class, getWorkflowPlan(),
                    htsfSample, false);
            clusterJob.setSiteName(siteName);
            clusterJob.addArgument(ClusterCLI.CLUSTERDIRECTORY, clusterDir.getAbsolutePath());
            graph.addVertex(clusterJob);
            graph.addEdge(parseClusterJob, clusterJob);

            // new job
            CondorJob mapspliceMultiThreadFusionJob = WorkflowJobFactory.createJob(++count,
                    MapSpliceMultiThreadCLI.class, getWorkflowPlan(), htsfSample);
            mapspliceMultiThreadFusionJob.setSiteName(siteName);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.FASTQ1,
                    remapUnmapped1RegexCatOutput.getAbsolutePath());
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.FASTQ2,
                    remapUnmapped2RegexCatOutput.getAbsolutePath());
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.THREADS, 8);
            mapspliceMultiThreadFusionJob.setNumberOfProcessors(8);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.QUALITYSCALE, QualityScaleType.phred33);
            File normalSAMMapspliceMultiThreadFusionOutput = new File(fusionDir, "normal.sam");
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.MAPSPLICEOUT,
                    normalSAMMapspliceMultiThreadFusionOutput.getAbsolutePath());
            unmappedOutput = new File(fusionDir, "fusion_original_unmapped");
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.OUTPUTUNMAPPED,
                    unmappedOutput.getAbsolutePath());
            mapspliceMultiThreadJobOutput = new File(fusionDir, "bowtie_original.output");
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.OUTPUT,
                    mapspliceMultiThreadJobOutput.getAbsolutePath());
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.OPTIMIZEREPEATS);
            File clusterRegion = new File(resultDir, "cluster.txt");
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.CLUSTERREGION,
                    clusterRegion.getAbsolutePath());
            File fusion = new File(fusionDir, "fusion_alignments_original.sam");
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.FUSION, fusion.getAbsolutePath());
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.JUNCTIONINDEX,
                    synIndexPrefix.getAbsolutePath());
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.MINIMUMINTRON, 50);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMINTRONSINGLE, 200000);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMINTRONDOUBLE, 200000);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.MINIMUMLENGTH, 25);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.SEGMENTLENGTH, 25);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.FASTQREADFORMAT);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMINSERTIONS, 6);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMDELETIONS, 6);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMAPPENDMISMATCHES, 3);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.REFERENCESEQUENCEDIRECTORY,
                    chromosomeDirectory);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.REPORTUPTONALIGNMENTSPERREAD, 40);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.SUPPRESSALIGNMENTSOVER, 40);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.SPLICEMISMATCHES, 1);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.CHROMOSOMESIZE,
                    chromosomeSizesOutput.getAbsolutePath());
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.INDEX, bowtieIndexDirectory);
            graph.addVertex(mapspliceMultiThreadFusionJob);
            graph.addEdge(clusterJob, mapspliceMultiThreadFusionJob);
            graph.addEdge(remapUnmapped1RegexCatJob, mapspliceMultiThreadFusionJob);
            graph.addEdge(remapUnmapped2RegexCatJob, mapspliceMultiThreadFusionJob);

            // new job
            CondorJob fusionSAM2JunctionFilterAnchorNewFormatJob = WorkflowJobFactory.createJob(++count,
                    FusionSAM2JunctionFilterAnchorNewFormatCLI.class, getWorkflowPlan(), htsfSample);
            fusionSAM2JunctionFilterAnchorNewFormatJob.setSiteName(siteName);
            File fusionSAM2JunctionFilterAnchorNewFormatOutput = new File(fusionDir, "original_fusion_junction.txt");
            fusionSAM2JunctionFilterAnchorNewFormatJob.addArgument(FusionSAM2JunctionFilterAnchorNewFormatCLI.JUNCTION,
                    fusionSAM2JunctionFilterAnchorNewFormatOutput.getAbsolutePath());
            fusionSAM2JunctionFilterAnchorNewFormatJob.addArgument(FusionSAM2JunctionFilterAnchorNewFormatCLI.SAM,
                    fusion.getAbsolutePath());
            fusionSAM2JunctionFilterAnchorNewFormatJob.addArgument(
                    FusionSAM2JunctionFilterAnchorNewFormatCLI.MINIMUMANCHOR, 1);
            fusionSAM2JunctionFilterAnchorNewFormatJob.addArgument(
                    FusionSAM2JunctionFilterAnchorNewFormatCLI.READLENGTH, determineReadLengthOutput.getAbsolutePath());
            fusionSAM2JunctionFilterAnchorNewFormatJob.addArgument(
                    FusionSAM2JunctionFilterAnchorNewFormatCLI.REFERENCESEQUENCEDIRECTORY, chromosomeDirectory);
            graph.addVertex(fusionSAM2JunctionFilterAnchorNewFormatJob);
            graph.addEdge(mapspliceMultiThreadFusionJob, fusionSAM2JunctionFilterAnchorNewFormatJob);

            // new job
            CondorJob filterOriginalFusionJob = WorkflowJobFactory.createJob(++count, FilterOriginalFusionCLI.class,
                    getWorkflowPlan(), htsfSample);
            filterOriginalFusionJob.setSiteName(siteName);
            filterOriginalFusionJob.addArgument(FilterOriginalFusionCLI.JUNCTION,
                    fusionSAM2JunctionFilterAnchorNewFormatOutput.getAbsolutePath());
            filterOriginalFusionJob.addArgument(FilterOriginalFusionCLI.MINIMUMLPQ, 0.3);
            filterOriginalFusionJob.addArgument(FilterOriginalFusionCLI.MINIMUMMISMATCHES, 2);
            File filteredJunctions = new File(fusionDir, "original_fusion_junction.filtered.txt");
            filterOriginalFusionJob.addArgument(FilterOriginalFusionCLI.FILTEREDJUNCTIONS,
                    filteredJunctions.getAbsolutePath());
            File remainingJunctions = new File(fusionDir, "original_fusion_junction.remained.txt");
            filterOriginalFusionJob.addArgument(FilterOriginalFusionCLI.REMAININGJUNCTIONS,
                    remainingJunctions.getAbsolutePath());
            graph.addVertex(filterOriginalFusionJob);
            graph.addEdge(fusionSAM2JunctionFilterAnchorNewFormatJob, filterOriginalFusionJob);

            // new job
            CondorJob junctionDBFusionJob = WorkflowJobFactory.createJob(++count, JunctionDBFusionCLI.class,
                    getWorkflowPlan(), htsfSample);
            junctionDBFusionJob.setSiteName(siteName);
            junctionDBFusionJob.addArgument(JunctionDBFusionCLI.FUSIONJUNCTION, remainingJunctions.getAbsolutePath());
            junctionDBFusionJob.addArgument(JunctionDBFusionCLI.JUNCTION, bestJunctionOutput.getAbsolutePath());
            junctionDBFusionJob.addArgument(JunctionDBFusionCLI.MAXIMUMANCHOR, 38);
            junctionDBFusionJob.addArgument(JunctionDBFusionCLI.MAXIMUMTHRESHOLDEACH, 20);
            junctionDBFusionJob.addArgument(JunctionDBFusionCLI.MAXIMUMTHRESHOLDTOTAL, 50);
            junctionDBFusionJob.addArgument(JunctionDBFusionCLI.MINIMUMANCHORWIDTH, 2);
            File junctionDBFusionOutput = new File(fusionDir, "fusion_synthetic_sequence.txt");
            junctionDBFusionJob.addArgument(JunctionDBFusionCLI.OUTPUT, junctionDBFusionOutput.getAbsolutePath());
            junctionDBFusionJob.addArgument(JunctionDBFusionCLI.REFERENCESEQUENCEDIRECTORY, chromosomeDirectory);
            graph.addVertex(junctionDBFusionJob);
            graph.addEdge(filterOriginalFusionJob, junctionDBFusionJob);

            // new job
            CondorJob synFusionIndexBowtieBuildJob = WorkflowJobFactory.createJob(++count, BowtieBuildCLI.class,
                    getWorkflowPlan(), htsfSample, false);
            synFusionIndexBowtieBuildJob.setSiteName(siteName);
            synFusionIndexBowtieBuildJob.addArgument(BowtieBuildCLI.INPUT, junctionDBFusionOutput.getAbsolutePath());
            File synFusionIndexPrefix = new File(fusionDir, "syn_fusion_idx_prefix");
            synFusionIndexBowtieBuildJob.addArgument(BowtieBuildCLI.PREFIX, synFusionIndexPrefix.getAbsolutePath());
            graph.addVertex(synFusionIndexBowtieBuildJob);
            graph.addEdge(junctionDBFusionJob, synFusionIndexBowtieBuildJob);

            // new job
            mapspliceMultiThreadFusionJob = WorkflowJobFactory.createJob(++count, MapSpliceMultiThreadCLI.class,
                    getWorkflowPlan(), htsfSample);
            mapspliceMultiThreadFusionJob.setSiteName(siteName);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.FASTQ1,
                    remapUnmapped1RegexCatOutput.getAbsolutePath());
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.FASTQ2,
                    remapUnmapped2RegexCatOutput.getAbsolutePath());
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.THREADS, 8);
            mapspliceMultiThreadFusionJob.setNumberOfProcessors(8);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.QUALITYSCALE, QualityScaleType.phred33);
            normalSAMMapspliceMultiThreadFusionOutput = new File(fusionDir, "normal.sam");
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.MAPSPLICEOUT,
                    normalSAMMapspliceMultiThreadFusionOutput.getAbsolutePath());
            unmappedOutput = new File(fusionDir, "fusion_unmapped");
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.OUTPUTUNMAPPED,
                    unmappedOutput.getAbsolutePath());
            mapspliceMultiThreadJobOutput = new File(fusionDir, "bowtie.output");
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.OUTPUT,
                    mapspliceMultiThreadJobOutput.getAbsolutePath());
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.OPTIMIZEREPEATS);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.MINIMUMINTRON, 50);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMINTRONSINGLE, 200000);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMINTRONDOUBLE, 200000);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.MINIMUMLENGTH, 25);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.SEGMENTLENGTH, 25);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.FASTQREADFORMAT);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMINSERTIONS, 6);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMDELETIONS, 6);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.MAXIMUMAPPENDMISMATCHES, 3);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.REFERENCESEQUENCEDIRECTORY,
                    chromosomeDirectory);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.REPORTUPTONALIGNMENTSPERREAD, 40);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.SUPPRESSALIGNMENTSOVER, 40);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.SPLICEMISMATCHES, 1);
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.CHROMOSOMESIZE,
                    chromosomeSizesOutput.getAbsolutePath());
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.CLUSTERREGION,
                    clusterRegion.getAbsolutePath());
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.JUNCTIONINDEX,
                    synIndexPrefix.getAbsolutePath());
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.FUSIONINDEX,
                    synFusionIndexPrefix.getAbsolutePath());
            fusion = new File(fusionDir, "fusion_alignments.sam");
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.FUSION, fusion.getAbsolutePath());
            mapspliceMultiThreadFusionJob.addArgument(MapSpliceMultiThreadCLI.INDEX, bowtieIndexDirectory);
            graph.addVertex(mapspliceMultiThreadFusionJob);
            graph.addEdge(synFusionIndexBowtieBuildJob, mapspliceMultiThreadFusionJob);

            // new job
            CondorJob sortJob = WorkflowJobFactory.createJob(++count, SortCLI.class, getWorkflowPlan(), htsfSample,
                    false);
            sortJob.setSiteName(siteName);
            File sortOutput = new File(fusionDir, "combined_fusion_normal.sam");
            sortJob.addArgument(SortCLI.OUTPUT, sortOutput.getAbsolutePath());
            sortJob.addArgument(SortCLI.BUFFERSIZE, 3500000);
            sortJob.addArgument(SortCLI.KEY, "1,1");
            sortJob.addArgument(SortCLI.TMPDIRECTORY, tmpDir.getAbsolutePath());
            sortJob.addArgument(SortCLI.INPUT, singleAlignments.getAbsolutePath());
            sortJob.addArgument(SortCLI.INPUT, fusionPairedAlignments.getAbsolutePath());
            sortJob.addArgument(SortCLI.INPUT, fusion.getAbsolutePath());
            graph.addVertex(sortJob);
            graph.addEdge(mapspliceMultiThreadFusionJob, sortJob);

            // new job
            CondorJob fusionUnmapped1RegexCatJob = WorkflowJobFactory.createJob(++count, RegexCatCLI.class,
                    getWorkflowPlan(), htsfSample, false);
            fusionUnmapped1RegexCatJob.setSiteName(siteName);
            fusionUnmapped1RegexCatJob.addArgument(RegexCatCLI.DIRECTORY, fusionDir.getAbsolutePath());
            File fusionUnmapped1RegexCatOutput = new File(fusionDir, "fusion_unmapped.1");
            fusionUnmapped1RegexCatJob.addArgument(RegexCatCLI.OUTPUT, fusionUnmapped1RegexCatOutput.getAbsolutePath());
            fusionUnmapped1RegexCatJob.addArgument(RegexCatCLI.REGEX, "^fusion_unmapped\\.[0-9]\\.1");
            graph.addVertex(fusionUnmapped1RegexCatJob);
            graph.addEdge(mapspliceMultiThreadFusionJob, fusionUnmapped1RegexCatJob);

            // new job
            CondorJob fusionUnmapped2RegexCatJob = WorkflowJobFactory.createJob(++count, RegexCatCLI.class,
                    getWorkflowPlan(), htsfSample, false);
            fusionUnmapped2RegexCatJob.setSiteName(siteName);
            fusionUnmapped2RegexCatJob.addArgument(RegexCatCLI.DIRECTORY, fusionDir.getAbsolutePath());
            File fusionUnmapped2RegexCatOutput = new File(fusionDir, "fusion_unmapped.2");
            fusionUnmapped2RegexCatJob.addArgument(RegexCatCLI.OUTPUT, fusionUnmapped2RegexCatOutput.getAbsolutePath());
            fusionUnmapped2RegexCatJob.addArgument(RegexCatCLI.REGEX, "^fusion_unmapped\\.[0-9]\\.2");
            graph.addVertex(fusionUnmapped2RegexCatJob);
            graph.addEdge(mapspliceMultiThreadFusionJob, fusionUnmapped2RegexCatJob);

            // new job
            CondorJob fusionUnmapped1ReadsToUnmappedSAMJob = WorkflowJobFactory.createJob(++count,
                    ReadsToUnmappedSAMCLI.class, getWorkflowPlan(), htsfSample, false);
            fusionUnmapped1ReadsToUnmappedSAMJob.setSiteName(siteName);
            fusionUnmapped1ReadsToUnmappedSAMJob.addArgument(ReadsToUnmappedSAMCLI.INPUT,
                    fusionUnmapped1RegexCatOutput.getAbsolutePath());
            File fusionUnmapped1ReadsToUnmappedSAMOutput = new File(fusionDir, "fusion_unmapped.1.sam");
            fusionUnmapped1ReadsToUnmappedSAMJob.addArgument(ReadsToUnmappedSAMCLI.OUTPUT,
                    fusionUnmapped1ReadsToUnmappedSAMOutput.getAbsolutePath());
            graph.addVertex(fusionUnmapped1ReadsToUnmappedSAMJob);
            graph.addEdge(fusionUnmapped1RegexCatJob, fusionUnmapped1ReadsToUnmappedSAMJob);

            // new job
            CondorJob fusionUnmapped2ReadsToUnmappedSAMJob = WorkflowJobFactory.createJob(++count,
                    ReadsToUnmappedSAMCLI.class, getWorkflowPlan(), htsfSample, false);
            fusionUnmapped2ReadsToUnmappedSAMJob.setSiteName(siteName);
            fusionUnmapped2ReadsToUnmappedSAMJob.addArgument(ReadsToUnmappedSAMCLI.INPUT,
                    fusionUnmapped2RegexCatOutput.getAbsolutePath());
            File fusionUnmapped2ReadsToUnmappedSAMOutput = new File(fusionDir, "fusion_unmapped.2.sam");
            fusionUnmapped2ReadsToUnmappedSAMJob.addArgument(ReadsToUnmappedSAMCLI.OUTPUT,
                    fusionUnmapped2ReadsToUnmappedSAMOutput.getAbsolutePath());
            graph.addVertex(fusionUnmapped2ReadsToUnmappedSAMJob);
            graph.addEdge(fusionUnmapped2RegexCatJob, fusionUnmapped2ReadsToUnmappedSAMJob);

            // new job
            alignmentHandlerMultiJob = WorkflowJobFactory.createJob(++count, AlignmentHandlerMultiCLI.class,
                    getWorkflowPlan(), htsfSample, false);
            alignmentHandlerMultiJob.setSiteName(siteName);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.ADDSOFTCLIP, 1);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.AVERAGEFRAGMENTLENGTH, 225);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.BOUNDARY, 36);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.CHROMOSOMESIZEFILE,
                    chromosomeSizesOutput.getAbsolutePath());
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.ENCOMPASSINGFUSIONREGIONEXTENSION, 50000);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.FILTEREDALIGNMENTAPPEND,
                    "_filtered_fusion_alignments");
            filteredAlignmentBase = new File(fusionDir, "_filtered_fusion_alignments");
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.FILTEREDALIGNMENTBASE,
                    filteredAlignmentBase.getAbsolutePath());
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.FILTERFLAG, 12 + 32 + 128 + 1024);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.FRAGMENTLENGTH, 400);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.FRAGMENTLENGTHSD, 100);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.INPUT, sortOutput.getAbsolutePath());
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.INTRONDISTANCESD, 500);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MATEDISTANCESD, 100);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MAXIMUMANCHORDIFF, 50);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MAXIMUMDELETION, 6);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MAXIMUMHITS, 40);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MAXIMUMMATEDISTANCE, 50000);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.READLENGTHPROPERTIES,
                    determineReadLengthOutput.getAbsolutePath());
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MINIMUMANCHOR, 0);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MINIMUMCOVERAGE, 0);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MINIMUMENCOMPASSCOUNT, 1);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MINIMUMENTROPY, -0.0001);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MINIMUMINSERTION, 6);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MINIMUMJUNCTIONANCHOR, 10);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.MINIMUMMISMATCH, 5);
            alignmentHandlerMultiJob.addArgument(AlignmentHandlerMultiCLI.THREADS, 8);
            alignmentHandlerMultiJob.setNumberOfProcessors(8);
            graph.addVertex(alignmentHandlerMultiJob);
            graph.addEdge(sortJob, alignmentHandlerMultiJob);
            graph.addEdge(determineReadLengthJob, alignmentHandlerMultiJob);

            // new job
            sortJob = WorkflowJobFactory.createJob(++count, SortCLI.class, getWorkflowPlan(), htsfSample, false);
            sortJob.setSiteName(siteName);
            sortOutput = new File(tmpDir, "combined_unmapped.sam");
            sortJob.addArgument(SortCLI.OUTPUT, sortOutput.getAbsolutePath());
            sortJob.addArgument(SortCLI.BUFFERSIZE, 3500000);
            sortJob.addArgument(SortCLI.KEY, "1,1");
            sortJob.addArgument(SortCLI.TMPDIRECTORY, tmpDir.getAbsolutePath());
            sortJob.addArgument(SortCLI.INPUT,
                    new File(remapDir, "_filtered_normal_alignments.unmapped").getAbsolutePath());
            sortJob.addArgument(SortCLI.INPUT, fusionUnmapped1ReadsToUnmappedSAMOutput.getAbsolutePath());
            sortJob.addArgument(SortCLI.INPUT, fusionUnmapped2ReadsToUnmappedSAMOutput.getAbsolutePath());
            sortJob.addArgument(SortCLI.INPUT,
                    new File(fusionDir, "_filtered_fusion_alignments.unmapped").getAbsolutePath());
            graph.addVertex(sortJob);
            graph.addEdge(alignmentHandlerMultiJob, sortJob);
            graph.addEdge(fusionUnmapped1ReadsToUnmappedSAMJob, sortJob);
            graph.addEdge(fusionUnmapped2ReadsToUnmappedSAMJob, sortJob);

            // new job
            CondorJob setUnmappedBitFlagJob = WorkflowJobFactory.createJob(++count, SetUnmappedBitFlagCLI.class,
                    getWorkflowPlan(), htsfSample, false);
            setUnmappedBitFlagJob.setSiteName(siteName);
            setUnmappedBitFlagJob.addArgument(SetUnmappedBitFlagCLI.UNMAPPEDSAM, sortOutput.getAbsolutePath());
            File unmappedBitFlag = new File(tmpDir, "combined_unmapped_setbitflag.sam");
            setUnmappedBitFlagJob.addArgument(SetUnmappedBitFlagCLI.UNMAPPEDSETBIT, unmappedBitFlag.getAbsolutePath());
            graph.addVertex(setUnmappedBitFlagJob);
            graph.addEdge(sortJob, setUnmappedBitFlagJob);

            // new job
            CondorJob finalAlignmentsHeadedCatJob = WorkflowJobFactory.createJob(++count, CatCLI.class,
                    getWorkflowPlan(), htsfSample, false);
            finalAlignmentsHeadedCatJob.setSiteName(siteName);
            File catOutput = new File(tmpDir, "final_alignments_headed.sam");
            finalAlignmentsHeadedCatJob.addArgument(CatCLI.OUTPUT, catOutput.getAbsolutePath());
            finalAlignmentsHeadedCatJob.addArgument(CatCLI.FILES, chromosomeHeadOutput.getAbsolutePath());
            finalAlignmentsHeadedCatJob.addArgument(CatCLI.FILES, pairedAlignments.getAbsolutePath());
            finalAlignmentsHeadedCatJob.addArgument(CatCLI.FILES, singleAlignments.getAbsolutePath());
            finalAlignmentsHeadedCatJob.addArgument(CatCLI.FILES, fusionPairedAlignments.getAbsolutePath());
            finalAlignmentsHeadedCatJob.addArgument(CatCLI.FILES, new File(remapDir,
                    "_filtered_normal_alignments.unmapped").getAbsolutePath());
            finalAlignmentsHeadedCatJob.addArgument(CatCLI.FILES,
                    fusionUnmapped1ReadsToUnmappedSAMOutput.getAbsolutePath());
            finalAlignmentsHeadedCatJob.addArgument(CatCLI.FILES,
                    fusionUnmapped2ReadsToUnmappedSAMOutput.getAbsolutePath());
            finalAlignmentsHeadedCatJob.addArgument(CatCLI.FILES, new File(fusionDir,
                    "_filtered_fusion_alignments.unmapped").getAbsolutePath());
            finalAlignmentsHeadedCatJob.addArgument(CatCLI.FILES, sortOutput.getAbsolutePath());
            graph.addVertex(finalAlignmentsHeadedCatJob);
            graph.addEdge(setUnmappedBitFlagJob, finalAlignmentsHeadedCatJob);

            // new job
            CondorJob samtoolsViewJob = WorkflowJobFactory.createJob(++count, SAMToolsViewCLI.class, getWorkflowPlan(),
                    htsfSample, false);
            samtoolsViewJob.setSiteName(siteName);
            File samtoolsViewOutput = new File(outputDirectory, fastqLaneRootName + ".bam");
            samtoolsViewJob.addArgument(SAMToolsViewCLI.OUTPUT, samtoolsViewOutput.getAbsolutePath());
            samtoolsViewJob.addArgument(SAMToolsViewCLI.SAMINPUTFORMAT);
            samtoolsViewJob.addArgument(SAMToolsViewCLI.BAMFORMAT);
            samtoolsViewJob.addArgument(SAMToolsViewCLI.INPUT, catOutput.getAbsolutePath());
            graph.addVertex(samtoolsViewJob);
            graph.addEdge(finalAlignmentsHeadedCatJob, samtoolsViewJob);

            // new job
            CondorJob picardAddOrReplaceReadGroupsJob = WorkflowJobFactory.createJob(++count,
                    PicardAddOrReplaceReadGroupsCLI.class, getWorkflowPlan(), htsfSample);
            picardAddOrReplaceReadGroupsJob.setSiteName(siteName);
            picardAddOrReplaceReadGroupsJob.addArgument(PicardAddOrReplaceReadGroupsCLI.INPUT,
                    samtoolsViewOutput.getAbsolutePath());
            File fixRGOutput = new File(outputDirectory, samtoolsViewOutput.getName().replace(".bam", ".fixed-rg.bam"));
            picardAddOrReplaceReadGroupsJob.addArgument(PicardAddOrReplaceReadGroupsCLI.OUTPUT,
                    fixRGOutput.getAbsolutePath());
            picardAddOrReplaceReadGroupsJob.addArgument(PicardAddOrReplaceReadGroupsCLI.SORTORDER,
                    PicardSortOrderType.COORDINATE.toString().toLowerCase());
            picardAddOrReplaceReadGroupsJob.addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPID, htsfSample.getId()
                    .toString());
            picardAddOrReplaceReadGroupsJob.addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPLIBRARY, sampleName);
            picardAddOrReplaceReadGroupsJob.addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPPLATFORM, sequencerRun
                    .getPlatform().getInstrument());
            picardAddOrReplaceReadGroupsJob.addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPPLATFORMUNIT,
                    sequencerRun.getPlatform().getInstrumentModel());
            picardAddOrReplaceReadGroupsJob
                    .addArgument(PicardAddOrReplaceReadGroupsCLI.READGROUPSAMPLENAME, sampleName);
            graph.addVertex(picardAddOrReplaceReadGroupsJob);
            graph.addEdge(samtoolsViewJob, picardAddOrReplaceReadGroupsJob);

            // new job
            CondorJob samtoolsSortJob = WorkflowJobFactory.createJob(++count, SAMToolsSortCLI.class, getWorkflowPlan(),
                    htsfSample);
            samtoolsSortJob.setSiteName(siteName);
            samtoolsSortJob.addArgument(SAMToolsSortCLI.INPUT, fixRGOutput.getAbsolutePath());
            File samtoolsSortOut = new File(outputDirectory, fixRGOutput.getName().replace(".bam", ".sorted.bam"));
            samtoolsSortJob.addArgument(SAMToolsSortCLI.OUTPUT, samtoolsSortOut.getAbsolutePath());
            graph.addVertex(samtoolsSortJob);
            graph.addEdge(picardAddOrReplaceReadGroupsJob, samtoolsSortJob);

            // new job
            CondorJob samtoolsIndexJob = WorkflowJobFactory.createJob(++count, SAMToolsIndexCLI.class,
                    getWorkflowPlan(), htsfSample);
            samtoolsIndexJob.setSiteName(siteName);
            samtoolsIndexJob.addArgument(SAMToolsIndexCLI.INPUT, samtoolsSortOut.getAbsolutePath());
            File samtoolsIndexOutput = new File(outputDirectory, samtoolsSortOut.getName().replace(".bam", ".bai"));
            samtoolsIndexJob.addArgument(SAMToolsIndexCLI.OUTPUT, samtoolsIndexOutput.getAbsolutePath());
            graph.addVertex(samtoolsIndexJob);
            graph.addEdge(samtoolsSortJob, samtoolsIndexJob);

            // new job
            CondorJob samtoolsFlagstatJob = WorkflowJobFactory.createJob(++count, SAMToolsFlagstatCLI.class,
                    getWorkflowPlan(), htsfSample);
            samtoolsFlagstatJob.setSiteName(siteName);
            samtoolsFlagstatJob.addArgument(SAMToolsFlagstatCLI.INPUT, samtoolsSortOut.getAbsolutePath());
            File samtoolsFlagstatOut = new File(outputDirectory, samtoolsSortOut.getName().replace(".bam", ".flagstat"));
            samtoolsFlagstatJob.addArgument(SAMToolsFlagstatCLI.OUTPUT, samtoolsFlagstatOut.getAbsolutePath());
            graph.addVertex(samtoolsFlagstatJob);
            graph.addEdge(samtoolsIndexJob, samtoolsFlagstatJob);

            // new job
            CondorJob ubuSamJunctionJob = WorkflowJobFactory.createJob(++count, UBUSamJunctionCLI.class,
                    getWorkflowPlan(), htsfSample);
            ubuSamJunctionJob.setSiteName(siteName);
            ubuSamJunctionJob.addArgument(UBUSamJunctionCLI.JUNCTIONS, junctions);
            ubuSamJunctionJob.addArgument(UBUSamJunctionCLI.INPUT, samtoolsSortOut.getAbsolutePath());
            File ubuSamJunctionOut = new File(outputDirectory, samtoolsSortOut.getName().replace(".bam",
                    ".junction_quantification.txt"));
            ubuSamJunctionJob.addArgument(UBUSamJunctionCLI.OUTPUT, ubuSamJunctionOut.getAbsolutePath());
            graph.addVertex(ubuSamJunctionJob);
            graph.addEdge(samtoolsIndexJob, ubuSamJunctionJob);

            // new job
            CondorJob coverageBedJob = WorkflowJobFactory.createJob(++count, CoverageBedCLI.class, getWorkflowPlan(),
                    htsfSample);
            coverageBedJob.setSiteName(siteName);
            coverageBedJob.addArgument(CoverageBedCLI.INPUT, samtoolsSortOut.getAbsolutePath());
            coverageBedJob.addArgument(CoverageBedCLI.BED, compositeExons);
            coverageBedJob.addArgument(CoverageBedCLI.SPLITBED);
            File coverageBedOut = new File(outputDirectory, samtoolsSortOut.getName().replace(".bam",
                    ".coverageBedOut.txt"));
            coverageBedJob.addArgument(CoverageBedCLI.OUTPUT, coverageBedOut.getAbsolutePath());
            graph.addVertex(coverageBedJob);
            graph.addEdge(samtoolsSortJob, coverageBedJob);

            // new job
            CondorJob normBedExonQuantJob = WorkflowJobFactory.createJob(++count, NormBedExonQuantCLI.class,
                    getWorkflowPlan(), htsfSample);
            normBedExonQuantJob.setSiteName(siteName);
            normBedExonQuantJob.addArgument(NormBedExonQuantCLI.INFILE, coverageBedOut.getAbsolutePath());
            normBedExonQuantJob.addArgument(NormBedExonQuantCLI.COMPOSITEBED, compositeExons);
            File normBedExonQuantOut = new File(outputDirectory, coverageBedOut.getName().replace(
                    ".coverageBedOut.txt", ".normBedExonQuantOut.txt"));
            normBedExonQuantJob.addArgument(NormBedExonQuantCLI.OUTFILE, normBedExonQuantOut.getAbsolutePath());
            graph.addVertex(normBedExonQuantJob);
            graph.addEdge(coverageBedJob, normBedExonQuantJob);

            // new job
            CondorJob sortBAMByReferenceAndNameJob = WorkflowJobFactory.createJob(++count,
                    SortByReferenceAndNameCLI.class, getWorkflowPlan(), htsfSample);
            sortBAMByReferenceAndNameJob.setSiteName(siteName);
            sortBAMByReferenceAndNameJob
                    .addArgument(SortByReferenceAndNameCLI.INPUT, samtoolsSortOut.getAbsolutePath());
            File sortBAMByReferenceAndNameOut = new File(outputDirectory, samtoolsSortOut.getName().replace(".bam",
                    ".refAndNameSort.bam"));
            sortBAMByReferenceAndNameJob.addArgument(SortByReferenceAndNameCLI.OUTPUT,
                    sortBAMByReferenceAndNameOut.getAbsolutePath());
            graph.addVertex(sortBAMByReferenceAndNameJob);
            graph.addEdge(samtoolsIndexJob, sortBAMByReferenceAndNameJob);

            // new job
            CondorJob ubuSamTranslateJob = WorkflowJobFactory.createJob(++count, UBUSamTranslateCLI.class,
                    getWorkflowPlan(), htsfSample);
            ubuSamTranslateJob.setSiteName(siteName);
            ubuSamTranslateJob.addArgument(UBUSamTranslateCLI.BED, bed);
            ubuSamTranslateJob.addArgument(UBUSamTranslateCLI.INPUT, sortBAMByReferenceAndNameOut.getAbsolutePath());
            File ubuSamTranslateOut = new File(outputDirectory, sortBAMByReferenceAndNameOut.getName().replace(".bam",
                    ".transcriptomeAlignments.bam"));
            ubuSamTranslateJob.addArgument(UBUSamTranslateCLI.OUTPUT, ubuSamTranslateOut.getAbsolutePath());
            ubuSamTranslateJob.addArgument(UBUSamTranslateCLI.ORDER,
                    getWorkflowBeanService().getAttributes().get("order"));
            ubuSamTranslateJob.addArgument(UBUSamTranslateCLI.XGTAGS);
            ubuSamTranslateJob.addArgument(UBUSamTranslateCLI.REVERSE);
            graph.addVertex(ubuSamTranslateJob);
            graph.addEdge(sortBAMByReferenceAndNameJob, ubuSamTranslateJob);

            // new job
            CondorJob ubuSamFilterJob = WorkflowJobFactory.createJob(++count, UBUSamFilterCLI.class, getWorkflowPlan(),
                    htsfSample);
            ubuSamFilterJob.setSiteName(siteName);
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
            CondorJob rsemJob = WorkflowJobFactory.createJob(++count, RSEMCalculateExpressionCLI.class,
                    getWorkflowPlan(), htsfSample);
            rsemJob.setSiteName(siteName);
            rsemJob.addArgument(RSEMCalculateExpressionCLI.PAIREDEND);
            rsemJob.addArgument(RSEMCalculateExpressionCLI.BAM);
            rsemJob.addArgument(RSEMCalculateExpressionCLI.ESTIMATERSPD);
            rsemJob.addArgument(RSEMCalculateExpressionCLI.THREADS, "4");
            rsemJob.setNumberOfProcessors(4);
            rsemJob.addArgument(RSEMCalculateExpressionCLI.BAMFILE, ubuSamFilterOut.getAbsolutePath());
            rsemJob.addArgument(RSEMCalculateExpressionCLI.REFERENCESEQUENCE, referenceSequencePrefix);
            File outputPrefix = new File(outputDirectory, "rsem");
            rsemJob.addArgument(RSEMCalculateExpressionCLI.OUTPUT, outputPrefix.getAbsolutePath());
            graph.addVertex(rsemJob);
            graph.addEdge(ubuSamFilterJob, rsemJob);

            // new job
            CondorJob rsemISOFormsResultStripTabJob = WorkflowJobFactory.createJob(++count, StripTrailingTabsCLI.class,
                    getWorkflowPlan(), htsfSample);
            rsemISOFormsResultStripTabJob.setSiteName(siteName);
            File isoformResults = new File(outputDirectory, outputPrefix.getName() + ".isoforms.results");
            rsemISOFormsResultStripTabJob.addArgument(StripTrailingTabsCLI.INPUT, isoformResults.getAbsolutePath());
            File isoformResultsStripped = new File(outputDirectory, outputPrefix.getName()
                    + ".isoforms.stripped.results");
            rsemISOFormsResultStripTabJob.addArgument(StripTrailingTabsCLI.OUTPUT,
                    isoformResultsStripped.getAbsolutePath());
            rsemISOFormsResultStripTabJob.setPostScript(String.format("/bin/mv %s %s",
                    isoformResultsStripped.getAbsolutePath(), isoformResults.getAbsolutePath()));
            graph.addVertex(rsemISOFormsResultStripTabJob);
            graph.addEdge(rsemJob, rsemISOFormsResultStripTabJob);

            // new job
            CondorJob pruneISOFormsFromGeneQuantFileJob = WorkflowJobFactory.createJob(++count,
                    PruneISOFormsFromGeneQuantFileCLI.class, getWorkflowPlan(), htsfSample);
            pruneISOFormsFromGeneQuantFileJob.setSiteName(siteName);
            File geneResultsFile = new File(outputDirectory, outputPrefix.getName() + ".genes.results");
            pruneISOFormsFromGeneQuantFileJob.addArgument(PruneISOFormsFromGeneQuantFileCLI.GENERESULTS,
                    geneResultsFile.getAbsolutePath());
            File origGeneResultsFile = new File(outputDirectory, "orig." + outputPrefix.getName() + ".genes.results");
            pruneISOFormsFromGeneQuantFileJob.addArgument(PruneISOFormsFromGeneQuantFileCLI.ORIGGENERESULTS,
                    origGeneResultsFile.getAbsolutePath());
            graph.addVertex(pruneISOFormsFromGeneQuantFileJob);
            graph.addEdge(rsemISOFormsResultStripTabJob, pruneISOFormsFromGeneQuantFileJob);

            // new job
            CondorJob normalizeGeneQuantJob = WorkflowJobFactory.createJob(++count, NormalizeQuartileCLI.class,
                    getWorkflowPlan(), htsfSample);
            normalizeGeneQuantJob.setSiteName(siteName);
            normalizeGeneQuantJob.addArgument(NormalizeQuartileCLI.COLUMN, "2");
            normalizeGeneQuantJob.addArgument(NormalizeQuartileCLI.QUANTILE, "75");
            normalizeGeneQuantJob.addArgument(NormalizeQuartileCLI.TARGET, "1000");
            normalizeGeneQuantJob.addArgument(NormalizeQuartileCLI.INPUT, geneResultsFile.getAbsolutePath());
            File normalizeGeneQuantOut = new File(outputDirectory, outputPrefix.getName() + ".genes.normalized_results");
            normalizeGeneQuantJob.addArgument(NormalizeQuartileCLI.OUTPUT, normalizeGeneQuantOut.getAbsolutePath());
            graph.addVertex(normalizeGeneQuantJob);
            graph.addEdge(pruneISOFormsFromGeneQuantFileJob, normalizeGeneQuantJob);

            // new job
            CondorJob normalizeISOFormQuantJob = WorkflowJobFactory.createJob(++count, NormalizeQuartileCLI.class,
                    getWorkflowPlan(), htsfSample);
            normalizeISOFormQuantJob.setSiteName(siteName);
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

            // new job
            CondorJob removeJob = WorkflowJobFactory.createJob(++count, RemoveCLI.class, getWorkflowPlan(), htsfSample);
            removeJob.setSiteName(siteName);
            removeJob.addArgument(RemoveCLI.FILE, gunzippedFastqR1.getAbsolutePath());
            removeJob.addArgument(RemoveCLI.FILE, gunzippedFastqR2.getAbsolutePath());
            removeJob.addArgument(RemoveCLI.FILE, tmpDir.getAbsolutePath());
            removeJob.addArgument(RemoveCLI.FILE, fastqFormatterR1Out.getAbsolutePath());
            removeJob.addArgument(RemoveCLI.FILE, fastqFormatterR2Out.getAbsolutePath());
            removeJob.addArgument(RemoveCLI.FILE, samtoolsViewOutput.getAbsolutePath());
            removeJob.addArgument(RemoveCLI.FILE, fixRGOutput.getAbsolutePath());
            graph.addVertex(removeJob);
            graph.addEdge(normalizeISOFormQuantJob, removeJob);
            graph.addEdge(normalizeGeneQuantJob, removeJob);

        }

        return graph;
    }
}
